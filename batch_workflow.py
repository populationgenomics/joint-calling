#!/usr/bin/env python3

"""
Hail Batch workflow to perform joint calling, sample QC, and variant QC with VQSR and 
random forest methods on a WGS germline callset.

The workflow is parametrised by the access level, the dataset name, 
batch names and the output version.

It must be only run with the CPG analysis-runner:
https://github.com/populationgenomics/analysis-runner (see helper script `driver_for_analysis_runner.sh` for analysis-runner submissions)
"""

import os
import subprocess
import tempfile
import traceback
from os.path import join, dirname, abspath, basename
from typing import List, Optional, Tuple, Collection
import logging
import click
import pandas as pd
import hailtop.batch as hb
from hailtop.batch.job import Job
from analysis_runner import dataproc
from sample_metadata import AnalysisApi, AnalysisModel

from joint_calling import utils
from joint_calling.utils import can_reuse
from joint_calling.variant_qc import add_variant_qc_jobs
from joint_calling import sm_utils

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
# @click.option('--dataset-version', 'dataset_version', type=str, required=True)
@click.option(
    '-n',
    '--namespace',
    'output_namespace',
    type=click.Choice(['main', 'test', 'tmp']),
    help='The bucket namespace to write the results to',
)
@click.option(
    '--analysis-project',
    'analysis_project',
    default='tob-wgs',
    help='Overrides the SM project name to write analysis entries to',
)
@click.option(
    '--input-project',
    'input_projects',
    multiple=True,
    help='Only read samples that belong to these project(s). Can be multiple.',
)
@click.option('--output-version', 'output_version', type=str, required=True)
@click.option(
    '--output-project',
    'output_projects',
    multiple=True,
    help='Only create reports for the project(s). Can be multiple. '
    'The name will be suffixed with the dataset version (set by --version)',
)
@click.option(
    '--use-gnarly-genotyper',
    'use_gnarly_genotyper',
    is_flag=True,
    help='Whether to use a GenomicsDB and Gnarly genotyper to merge GVCFs instead of '
    'the Hail combiner',
)
@click.option(
    '--ped-file',
    'ped_file',
    help='PED file with family information',
    type=str,
    callback=utils.get_validation_callback(ext='ped', must_exist=True),
)
@click.option(
    '--filter-cutoffs-file',
    'filter_cutoffs_path',
    callback=utils.get_validation_callback(ext='yaml', must_exist=True),
    help=f'YAML file with filtering cutoffs. '
    f'Default is the file within the package: {utils.get_filter_cutoffs()}',
)
@click.option('--keep-scratch', 'keep_scratch', is_flag=True)
@click.option(
    '--reuse-scratch-run-id',
    'reuse_scratch_run_id',
    help='Run ID to reuse scratch from',
)
@click.option(
    '--overwrite/--reuse',
    'overwrite',
    is_flag=True,
    help='if an intermediate or a final file exists, skip running the code '
    'that generates it.',
)
@click.option(
    '--skip-sample',
    '-S',
    'skip_samples',
    multiple=True,
    help='Don\'t process specified samples. Can be set multiple times.',
)
@click.option(
    '--scatter-count',
    'scatter_count',
    type=int,
    required=True,
    help='Number of secondary workers for Dataproc clusters, as well as the'
    'number of shards to parition data for the AS-VQSR analysis',
)
@click.option(
    '--pca_pop',
    'pca_pop',
    help='if specified, a separate PCA will be produced with samples trained only'
    'on this population',
)
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option('--run-vqsr/--skip-vqsr', 'run_vqsr', is_flag=True, default=True)
@click.option('--run-rf/--skip-rf', 'run_rf', is_flag=True, default=False)
def main(  # pylint: disable=too-many-arguments,too-many-locals,too-many-statements
    output_namespace: str,
    analysis_project: str,
    input_projects: List[str],  # pylint: disable=unused-argument
    output_version: str,
    output_projects: Optional[List[str]],  # pylint: disable=unused-argument
    use_gnarly_genotyper: bool,
    ped_file: str,
    filter_cutoffs_path: str,
    keep_scratch: bool,
    reuse_scratch_run_id: str,  # pylint: disable=unused-argument
    overwrite: bool,
    skip_samples: Collection[str],
    scatter_count: int,
    pca_pop: Optional[str],
    dry_run: bool,
    run_vqsr: bool,
    run_rf: bool,
):  # pylint: disable=missing-function-docstring
    # Determine bucket paths
    if output_namespace in ['test', 'tmp']:
        tmp_bucket_suffix = 'test-tmp'
    else:
        tmp_bucket_suffix = 'main-tmp'

    if output_namespace in ['test', 'main']:
        output_suffix = output_namespace
        output_analysis_suffix = f'{output_namespace}-analysis'
        web_bucket_suffix = f'{output_namespace}-web'
    else:
        output_suffix = 'test-tmp'
        output_analysis_suffix = 'test-tmp'
        web_bucket_suffix = 'test-tmp'

    ptrn = f'gs://cpg-{analysis_project}-{{suffix}}/joint-calling/{output_version}'
    work_bucket = ptrn.format(suffix=tmp_bucket_suffix)
    output_analysis_bucket = ptrn.format(suffix=output_analysis_suffix)
    web_bucket = ptrn.format(suffix=web_bucket_suffix)

    combiner_bucket = f'{work_bucket}/combiner'
    sample_qc_bucket = f'{work_bucket}/sample_qc'

    mt_output_bucket = f'gs://cpg-{analysis_project}-{output_suffix}/mt'
    raw_combined_mt_path = f'{combiner_bucket}/{output_version}-raw.mt'
    filtered_combined_mt_path = f'{mt_output_bucket}/{output_version}.mt'
    filtered_combined_nonref_mt_path = f'{mt_output_bucket}/{output_version}-nonref.mt'
    _create_mt_readme(
        raw_combined_mt_path,
        filtered_combined_mt_path,
        filtered_combined_nonref_mt_path,
        mt_output_bucket,
    )

    # Initialize Hail Batch
    billing_project = os.getenv('HAIL_BILLING_PROJECT') or analysis_project
    hail_bucket = os.environ.get('HAIL_BUCKET')
    if not hail_bucket or keep_scratch or reuse_scratch_run_id:
        # Scratch files are large, so we want to use the temporary bucket for them
        hail_bucket = f'{work_bucket}/hail'
    logger.info(
        f'Starting hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=hail_bucket.replace('gs://', ''),
    )
    b = hb.Batch(
        f'Joint calling: {analysis_project}'
        f', version: {output_version}'
        f', projects: {output_projects}'
        + (f', limit input projects to {output_projects}' if output_projects else '')
        + f', namespace: {output_namespace}',
        backend=backend,
    )

    # Scripts path to pass to dataproc when submitting scripts
    scripts_dir = abspath(join(dirname(__file__), 'scripts'))

    filter_cutoffs_d = utils.get_filter_cutoffs(filter_cutoffs_path)

    samples_df, samples_csv_path, pre_combiner_jobs = _add_pre_combiner_jobs(
        b=b,
        work_bucket=join(work_bucket, 'pre-combine'),
        output_bucket=combiner_bucket,
        overwrite=overwrite,
        input_projects=input_projects,
        analysis_project=analysis_project,
        is_test=output_namespace in ['test', 'tmp'],
        skip_samples=skip_samples,
    )

    if use_gnarly_genotyper:
        combiner_job = b.new_job('Not implemented')
        # combiner_job, vcf_path = make_joint_genotype_jobs(
        #     b=b,
        #     genomicsdb_bucket=join(work_bucket, 'genomicsdbs'),
        #     samples_df=samples_df,
        #     reference=reference,
        #     dbsnp=utils.DBSNP_VCF,
        #     out_bucket=out_bucket,
        #     tmp_bucket=out_bucket,
        #     local_tmp_dir=local_tmp_dir,
        #     overwrite=overwrite,
        #     analysis_project=analysis_project,
        #     completed_analysis=jc_analysis,
        #     depends_on=gvcf_jobs,
        # )
    else:
        if not can_reuse(raw_combined_mt_path, overwrite):
            combiner_job = dataproc.hail_dataproc_job(
                b,
                f'{scripts_dir}/combine_gvcfs.py '
                f'--meta-csv {samples_csv_path} '
                f'--out-mt {raw_combined_mt_path} '
                f'--bucket {combiner_bucket}/work '
                f'--hail-billing {billing_project} '
                f'--n-partitions {scatter_count * 25}',
                max_age='8h',
                packages=utils.DATAPROC_PACKAGES,
                num_secondary_workers=scatter_count,
                depends_on=pre_combiner_jobs,
                job_name='Combine GVCFs',
            )
        else:
            combiner_job = b.new_job('Combine GVCFs [reuse]')

    sample_qc_job, hard_filter_ht_path, meta_ht_path = _add_sample_qc_jobs(
        b=b,
        mt_path=raw_combined_mt_path,
        samples_csv_path=samples_csv_path,
        work_bucket=sample_qc_bucket,
        out_analysis_bucket=output_analysis_bucket,
        filter_cutoffs_path=filter_cutoffs_path,
        scripts_dir=scripts_dir,
        overwrite=overwrite,
        scatter_count=scatter_count,
        depends_on=[combiner_job],
        billing_project=billing_project,
        sample_count=len(samples_df),
        age_csv_path=f'gs://cpg-{analysis_project}-main-analysis/metadata/age.csv',
        is_test=output_namespace in ['test', 'tmp'],
        pca_pop=pca_pop,
    )

    if run_rf or run_vqsr:
        var_qc_job, vqsr_final_filter_ht = add_variant_qc_jobs(
            b=b,
            work_bucket=join(work_bucket, 'variant_qc'),
            web_bucket=join(web_bucket, 'variant_qc'),
            raw_combined_mt_path=raw_combined_mt_path,
            hard_filter_ht_path=hard_filter_ht_path,
            meta_ht_path=meta_ht_path,
            samples_df=samples_df,
            scripts_dir=scripts_dir,
            ped_file=ped_file,
            overwrite=overwrite,
            run_rf=run_rf,
            vqsr_params_d=filter_cutoffs_d['vqsr'],
            scatter_count=scatter_count,
            depends_on=[sample_qc_job],
        )
        if not can_reuse(filtered_combined_mt_path, overwrite):
            var_qc_job = dataproc.hail_dataproc_job(
                b,
                f'{scripts_dir}/make_finalised_mt.py --overwrite '
                f'--mt {raw_combined_mt_path} '
                f'--final-filter-ht {vqsr_final_filter_ht} '
                f'--out-mt {filtered_combined_mt_path} '
                f'--out-nonref-mt {filtered_combined_nonref_mt_path} '
                f'--meta-ht {meta_ht_path} ',
                max_age='8h',
                packages=utils.DATAPROC_PACKAGES,
                num_secondary_workers=scatter_count,
                depends_on=[var_qc_job, combiner_job],
                job_name='Making final MT',
            )
        else:
            var_qc_job = b.new_job('Making final MT [reuse]')
    else:
        var_qc_job = b.new_job('Var QC [skip]')

    # Interacting with the sample metadata server.
    aapi = AnalysisApi()
    # 1. Create a "queued" analysis
    am = AnalysisModel(
        type='joint-calling',
        output=filtered_combined_nonref_mt_path,
        status='queued',
        sample_ids=samples_df.s,
    )
    try:
        aid = aapi.create_new_analysis(project=analysis_project, analysis_model=am)
        # 2. Queue a job that updates the status to "in-progress"
        sm_in_progress_j = sm_utils.make_sm_in_progress_job(
            b, 'joint-calling', analysis_project, aid
        )
        # 2. Queue a job that updates the status to "completed"
        sm_completed_j = sm_utils.make_sm_completed_job(
            b, 'joint-calling', analysis_project, aid
        )
        # Set up dependencies
        combiner_job.depends_on(sm_in_progress_j)
        if pre_combiner_jobs:
            sm_in_progress_j.depends_on(*pre_combiner_jobs)
        sm_completed_j.depends_on(var_qc_job)
        logger.info(f'Queueing {am.type} with analysis ID: {aid}')
    except Exception:  # pylint: disable=broad-except
        print(traceback.format_exc())

    b.run(dry_run=dry_run, delete_scratch_on_exit=not keep_scratch, wait=False)


def _create_mt_readme(
    raw_combined_mt_path: str,
    filtered_combined_mt_path: str,
    filtered_combined_nonref_mt_path: str,
    mt_output_bucket: str,
):
    local_tmp_dir = tempfile.mkdtemp()
    readme_local_fpath = join(local_tmp_dir, 'README.txt')
    readme_gcs_fpath = join(mt_output_bucket, 'README.txt')
    with open(readme_local_fpath, 'w') as f:
        f.write(f'Unfiltered:\n{raw_combined_mt_path}')
        f.write(
            f'AS-VQSR soft-filtered\n(use `mt.filter_rows(hl.is_missing(mt.filters)` '
            f'to hard-filter):\n {basename(filtered_combined_mt_path)}'
        )
        f.write(
            f'AS-VQSR soft-filtered, without reference blocks:\n'
            f'{basename(filtered_combined_nonref_mt_path)}'
        )
    subprocess.run(['gsutil', 'cp', readme_local_fpath, readme_gcs_fpath], check=False)


def _add_sample_qc_jobs(
    b: hb.Batch,
    mt_path: str,
    samples_csv_path: str,
    work_bucket: str,
    out_analysis_bucket: str,
    age_csv_path: str,
    filter_cutoffs_path: Optional[str],
    scripts_dir: str,
    overwrite: bool,
    scatter_count: int,
    sample_count: int,  # pylint: disable=unused-argument
    depends_on: Optional[List[Job]] = None,
    billing_project: Optional[str] = None,
    is_test: bool = False,
    pca_pop: Optional[str] = None,
) -> Tuple[Job, str, str]:

    sample_qc_bucket = work_bucket

    pca_with_hgdp_mt_path = join(work_bucket, 'mt_subset_for_pca_with_hgdp.mt')
    pca_mt_path = join(work_bucket, 'mt_subset_for_pca.mt')
    if not all(can_reuse(fp, overwrite) for fp in [pca_mt_path, pca_with_hgdp_mt_path]):
        pcrelate_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/sample_qc_subset_mt_for_pca.py '
            + (f'--overwrite ' if overwrite else '')
            + f'--mt {mt_path} '
            f'--out-hgdp-union-mt {pca_with_hgdp_mt_path} '
            f'--out-mt {pca_mt_path} '
            + ('--is-test ' if is_test else '')
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=(depends_on or []),
            job_name='Subset MT for PCA',
        )
    else:
        pcrelate_job = b.new_job(f'Subset MT for PCA [reuse]')

    relatedness_ht_path = join(out_analysis_bucket, 'relatedness.ht')
    if not can_reuse(relatedness_ht_path, overwrite):
        pcrelate_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/sample_qc_pcrelate.py '
            + (f'--overwrite ' if overwrite else '')
            + f'--pca-mt {pca_mt_path} '
            f'--out-relatedness-ht {relatedness_ht_path} '
            f'--tmp-bucket {work_bucket}/tmp '
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            # Spark would have problems shuffling on pre-emptible workers
            # (throw SparkException: Job aborted due to stage failure: ShuffleMapStage)
            # so we use num_workers instead of num_secondary_workers here
            num_workers=scatter_count,
            depends_on=[pcrelate_job],
            job_name='Run pcrelate',
        )
    else:
        pcrelate_job = b.new_job(f'Run pcrelate [reuse]')

    if filter_cutoffs_path:
        gcs_path = join(work_bucket, 'filter-cutoffs.yaml')
        subprocess.run(['gsutil', 'cp', filter_cutoffs_path, gcs_path], check=False)
        filter_cutoffs_param = f'--filter-cutoffs-file {gcs_path}'
    else:
        filter_cutoffs_param = ''

    hard_filter_ht_path = join(sample_qc_bucket, 'hard_filters.ht')
    pop_ht_path = join(sample_qc_bucket, 'pop.ht')
    regressed_metrics_ht_path = join(sample_qc_bucket, 'regress_metrics.ht')
    if any(
        not can_reuse(fp, overwrite)
        for fp in [
            hard_filter_ht_path,
            pop_ht_path,
            regressed_metrics_ht_path,
        ]
    ):
        sample_qc_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/sample_qc_calculate.py '
            + (f'--overwrite ' if overwrite else '')
            + f'{filter_cutoffs_param} '
            + f'--mt {mt_path} '
            f'--pca-hgdp-mt-path {pca_with_hgdp_mt_path} '
            f'--meta-csv {samples_csv_path} '
            f'--relatedness-ht {relatedness_ht_path} '
            f'--out-hard-filtered-samples-ht {hard_filter_ht_path} '
            f'--out-bucket {work_bucket} '
            f'--out-analysis-bucket {out_analysis_bucket} '
            f'--tmp-bucket {work_bucket}/tmp '
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='12h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=(depends_on or []) + [pcrelate_job],
            job_name='Sample QC',
        )
    else:
        sample_qc_job = b.new_job(f'Sample QC [reuse]')

    meta_ht_path = join(out_analysis_bucket, 'meta.ht')
    meta_tsv_path = join(out_analysis_bucket, 'meta.tsv')
    if any(
        not can_reuse(fp, overwrite)
        for fp in [
            meta_ht_path,
            meta_tsv_path,
        ]
    ):
        if utils.file_exists(age_csv_path):
            age_csv_param = f'--age-csv {age_csv_path} '
        else:
            age_csv_param = ''
        metadata_qc_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/sample_qc_write_metadata.py '
            + (f'--overwrite ' if overwrite else '')
            + f'{filter_cutoffs_param} '
            + f'--qc-bucket {work_bucket} '
            f'--tmp-bucket {work_bucket}/tmp '
            f'{age_csv_param}'
            f'--hard-filtered-samples-ht {hard_filter_ht_path} '
            f'--relatedness-ht {relatedness_ht_path} '
            f'--out-meta-ht {meta_ht_path} '
            f'--out-meta-tsv {meta_tsv_path} '
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            depends_on=[sample_qc_job],
            job_name=f'Write sample QC metadata',
        )
    else:
        metadata_qc_job = b.new_job(f'Write sample QC metadata [reuse]')

    if pca_pop:
        pop_pca_scores_ht_path = join(
            out_analysis_bucket,
            utils.PCA_SCORES_HT_NAME.replace('.ht', f'_{pca_pop}.ht'),
        )
        if not can_reuse(pop_pca_scores_ht_path, overwrite):
            dataproc.hail_dataproc_job(
                b,
                f'{scripts_dir}/pca_pop.py '
                + f'--pca-hgdp-mt-path {pca_with_hgdp_mt_path} '
                f'--pop {pca_pop} '
                f'--out-analysis-bucket {out_analysis_bucket} '
                f'--meta-ht {meta_ht_path} '
                f'--tmp-bucket {work_bucket}/tmp '
                + (f'--overwrite ' if overwrite else '')
                + (f'--hail-billing {billing_project} ' if billing_project else ''),
                max_age='8h',
                packages=utils.DATAPROC_PACKAGES,
                depends_on=[metadata_qc_job],
                job_name=f'Pop-specific PCA',
            )
        else:
            b.new_job(f'Pop-specific PCA [reuse]')

    return metadata_qc_job, hard_filter_ht_path, meta_ht_path


def _add_pre_combiner_jobs(
    b: hb.Batch,
    work_bucket: str,
    output_bucket: str,
    overwrite: bool,
    input_projects: List[str],
    analysis_project: str,  # pylint: disable=unused-argument
    is_test: bool,
    skip_samples: Optional[Collection[str]] = None,
) -> Tuple[pd.DataFrame, str, List[Job]]:
    """
    Add jobs that prepare GVCFs for the combiner, if needed.

    :param work_bucket: bucket to write intermediate files to
    :param output_bucket: bucket to write the GVCF combiner inputs to
    :param overwrite: ignore existing intermediate files
    :param is_test: read gvcfs from a test bucket instead of main
    :return: a Tuple of: a pandas dataframe with the sample metadata, a CSV file
    corresponding to that dataframe, and a list of jobs to wait for before
    submitting the combiner job
    """

    # File with the pointers to GVCFs to process along with metdata.
    # If it doesn't exist, we trigger a utuils.find_inputs(combiner_bucket) function
    # to find the GVCFs and the metadata given the requested batch ids.
    input_samples_tsv_path = join(work_bucket, 'samples.tsv')
    # Raw GVCFs need pre-processing before passing to the combiner. If the following
    # file exists, we assume the samples are pre-processed; otherwise, we add Batch
    # jobs to do the pre-processing.
    combiner_ready_samples_tsv_path = join(output_bucket, 'samples.tsv')
    subset_gvcf_jobs: List[Job] = []
    if can_reuse(combiner_ready_samples_tsv_path, overwrite):
        logger.info(
            f'Reading existing combiner-read inputs TSV {combiner_ready_samples_tsv_path}'
        )
        samples_df = pd.read_csv(combiner_ready_samples_tsv_path, sep='\t').set_index(
            's', drop=False
        )
        samples_df = samples_df[pd.notnull(samples_df.s)]
    else:
        if can_reuse(input_samples_tsv_path, overwrite):
            logger.info(f'Reading existing inputs TSV {input_samples_tsv_path}')
            samples_df = pd.read_csv(input_samples_tsv_path, sep='\t').set_index(
                's', drop=False
            )
        else:
            logger.info(
                f'Querying samples from the sample-metadata server '
                f'for the projects: {", ".join(input_projects)}'
            )
            samples_df = sm_utils.find_inputs(
                input_projects,
                is_test=is_test,
                skip_samples=skip_samples,
            )
            samples_df.to_csv(
                input_samples_tsv_path, index=False, sep='\t', na_rep='NA'
            )

        samples_df = samples_df[pd.notnull(samples_df.s)]
        subset_gvcf_jobs, samples_df = _add_prep_gvcfs_for_combiner_steps(
            b=b,
            samples_df=samples_df,
            output_gvcf_bucket=join(output_bucket, 'gvcf'),
        )
        samples_df.to_csv(
            combiner_ready_samples_tsv_path, index=False, sep='\t', na_rep='NA'
        )
        logger.info(
            f'Saved metadata with updated GVCFs to '
            f'{combiner_ready_samples_tsv_path}'
        )

    return samples_df, combiner_ready_samples_tsv_path, subset_gvcf_jobs


def _add_prep_gvcfs_for_combiner_steps(
    b,
    samples_df: pd.DataFrame,
    output_gvcf_bucket: str,
    depends_on: Optional[List[Job]] = None,
) -> Tuple[List[Job], pd.DataFrame]:
    """
    Add steps required to prepare GVCFs from combining
    """
    noalt_regions = b.read_input(utils.NOALT_REGIONS)

    jobs = []
    logger.info(f'Samples DF: {samples_df}')
    for s_id, external_id, gvcf_path in zip(
        samples_df.s, samples_df.external_id, samples_df.gvcf
    ):
        logger.info(
            f'Adding reblock and subset jobs for sample {s_id}, gvcf {gvcf_path}'
        )
        gvcf = b.read_input_group(
            **{'g.vcf.gz': gvcf_path, 'g.vcf.gz.tbi': gvcf_path + '.tbi'}
        )
        reblock_j = _add_reblock_gvcfs_step(b, gvcf, depends_on=depends_on)
        output_gvcf_path = join(output_gvcf_bucket, f'{s_id}.g.vcf.gz')
        jobs.append(
            _add_subset_noalt_step(
                b,
                input_gvcf=reblock_j.output_gvcf,
                output_gvcf_path=output_gvcf_path,
                noalt_regions=noalt_regions,
                depends_on=depends_on,
                external_sample_id=external_id,
                internal_sample_id=s_id,
            )
        )
        assert s_id
        samples_df.loc[s_id, ['gvcf']] = output_gvcf_path
        logger.info(f'Updating sample {s_id} gvcf to {output_gvcf_path}')
    logger.info(f'Updated sample DF: {samples_df}')
    return jobs, samples_df


def _add_reblock_gvcfs_step(
    b: hb.Batch,
    input_gvcf: hb.ResourceGroup,
    depends_on: Optional[List[Job]] = None,
) -> Job:
    """
    Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration
    """
    j = b.new_job('ReblockGVCF')
    j.image(utils.GATK_IMAGE)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'30G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )
    if depends_on:
        j.depends_on(*depends_on)

    j.command(
        f"""
    gatk --java-options "-Xms{mem_gb - 1}g" \\
        ReblockGVCF \\
        -V {input_gvcf['g.vcf.gz']} \\
        -do-qual-approx \\
        -O {j.output_gvcf['g.vcf.gz']} \\
        --create-output-variant-index true"""
    )
    return j


def _add_subset_noalt_step(
    b: hb.Batch,
    input_gvcf: hb.ResourceGroup,
    output_gvcf_path: str,
    noalt_regions: hb.ResourceFile,
    external_sample_id: str,
    internal_sample_id: str,
    depends_on: Optional[List[Job]] = None,
) -> Job:
    """
    1. Subset GVCF to main chromosomes to avoid downstream errors
    2. Removes the DS INFO field that is added to some HGDP GVCFs to avoid errors
       from Hail about mismatched INFO annotations
    3. Renames sample name from external_sample_id to internal_sample_id
    """
    j = b.new_job('SubsetToNoalt')
    j.image('quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2')
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'30G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )
    if depends_on:
        j.depends_on(*depends_on)
    j.command(
        f"""set -e

    bcftools view {input_gvcf['g.vcf.gz']} -T {noalt_regions} \\
        | bcftools annotate -x INFO/DS \\
        | bcftools reheader -s <(echo "{external_sample_id} {internal_sample_id}") \\
        | bcftools view -Oz -o {j.output_gvcf['g.vcf.gz']}

    bcftools index --tbi {j.output_gvcf['g.vcf.gz']}
        """
    )
    b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
    return j


if __name__ == '__main__':
    main()  # pylint: disable=E1120
