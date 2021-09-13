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
from os.path import join, basename
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
from joint_calling import pre_combiner

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
    '--num-ancestry-pcs',
    'num_ancestry_pcs',
    type=int,
    default=20,
)
@click.option(
    '--pca-pop',
    'pca_pop',
    help='if specified, a separate PCA will be produced with samples trained only'
    'on this population',
)
@click.option(
    '--pre-computed-hgdp-union-mt',
    'pre_computed_hgdp_unuon_mt_path',
    help='Useful for tests to save time on subsetting the large gnomAD matrix table',
)
@click.option('--dry-run', 'dry_run', is_flag=True)
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
    num_ancestry_pcs: int,
    pre_computed_hgdp_unuon_mt_path: str,
    dry_run: bool,
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
        output_analysis_suffix = 'test-tmp/analysis'
        web_bucket_suffix = 'test-tmp/web'

    ptrn = f'gs://cpg-{analysis_project}-{{suffix}}/joint-calling/{output_version}'
    work_bucket = ptrn.format(suffix=tmp_bucket_suffix)
    out_analysis_bucket = ptrn.format(suffix=output_analysis_suffix)
    web_bucket = ptrn.format(suffix=web_bucket_suffix)

    pre_combiner_bucket = f'{work_bucket}/pre_combiner'
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

    samples_df = find_inputs(
        output_bucket=work_bucket,
        overwrite=overwrite,
        input_projects=input_projects,
        is_test=output_namespace in ['test', 'tmp'],
        skip_samples=skip_samples,
    )

    (
        samples_df,
        samples_csv_path,
        pre_combiner_jobs,
    ) = pre_combiner.add_pre_combiner_jobs(
        b=b,
        samples_df=samples_df,
        pre_combiner_bucket=pre_combiner_bucket,
        output_suffix=output_suffix,
        overwrite=overwrite,
        analysis_project=analysis_project,
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
                f'{utils.SCRIPTS_DIR}/combine_gvcfs.py '
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
        sample_qc_bucket=sample_qc_bucket,
        out_analysis_bucket=out_analysis_bucket,
        out_web_bucket=web_bucket,
        filter_cutoffs_path=filter_cutoffs_path,
        overwrite=overwrite,
        scatter_count=scatter_count,
        combiner_job=combiner_job,
        billing_project=billing_project,
        sample_count=len(samples_df),
        age_csv_path=f'gs://cpg-{analysis_project}-main-analysis/metadata/age.csv',
        num_ancestry_pcs=num_ancestry_pcs,
        pca_pop=pca_pop,
        is_test=output_namespace in ['test', 'tmp'],
        pre_computed_hgdp_unuon_mt_path=pre_computed_hgdp_unuon_mt_path,
    )

    var_qc_job, vqsr_final_filter_ht = add_variant_qc_jobs(
        b=b,
        work_bucket=join(work_bucket, 'variant_qc'),
        web_bucket=join(web_bucket, 'variant_qc'),
        raw_combined_mt_path=raw_combined_mt_path,
        hard_filter_ht_path=hard_filter_ht_path,
        meta_ht_path=meta_ht_path,
        samples_df=samples_df,
        ped_file=ped_file,
        overwrite=overwrite,
        vqsr_params_d=utils.get_filter_cutoffs(filter_cutoffs_path)['vqsr'],
        scatter_count=scatter_count,
        depends_on=[sample_qc_job],
    )

    if not can_reuse(filtered_combined_mt_path, overwrite):
        var_qc_job = dataproc.hail_dataproc_job(
            b,
            f'{utils.SCRIPTS_DIR}/make_finalised_mt.py --overwrite '
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
    sample_qc_bucket: str,
    out_analysis_bucket: str,
    out_web_bucket: str,
    age_csv_path: str,
    filter_cutoffs_path: Optional[str],
    overwrite: bool,
    scatter_count: int,
    sample_count: int,  # pylint: disable=unused-argument
    combiner_job: Job,
    num_ancestry_pcs: int,
    pca_pop: Optional[str] = None,
    billing_project: Optional[str] = None,
    is_test: bool = False,
    pre_computed_hgdp_unuon_mt_path: Optional[str] = None,
) -> Tuple[Job, str, str]:

    if filter_cutoffs_path:
        gcs_path = join(sample_qc_bucket, 'filter-cutoffs.yaml')
        subprocess.run(['gsutil', 'cp', filter_cutoffs_path, gcs_path], check=False)
        filter_cutoffs_param = f'--filter-cutoffs-file {gcs_path}'
    else:
        filter_cutoffs_param = ''

    input_metadata_ht_path = join(sample_qc_bucket, 'input_metadata.ht')
    hard_filtered_samples_ht_path = join(sample_qc_bucket, 'hard_filtered_samples.ht')
    sex_ht_path = join(sample_qc_bucket, 'sex.ht')
    hail_sample_qc_ht_path = join(sample_qc_bucket, 'hail_sample_qc.ht')
    custom_qc_ht_path = join(sample_qc_bucket, 'custom_qc.ht')
    if not can_reuse(
        [
            input_metadata_ht_path,
            hard_filtered_samples_ht_path,
            sex_ht_path,
            hail_sample_qc_ht_path,
            custom_qc_ht_path,
        ],
        overwrite,
    ):
        sample_qc_hardfilter_job = dataproc.hail_dataproc_job(
            b,
            f'{utils.SCRIPTS_DIR}/sample_qc_hard_filters.py '
            f'--mt {mt_path} '
            f'--meta-csv {samples_csv_path} '
            f'{filter_cutoffs_param} '
            f'--out-input-metadata-ht {input_metadata_ht_path} '
            f'--out-hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
            f'--out-sex-ht {sex_ht_path} '
            f'--out-hail-sample-qc-ht {hail_sample_qc_ht_path} '
            f'--out-custom-qc-ht {custom_qc_ht_path} '
            f'--tmp-bucket {sample_qc_bucket}/tmp '
            + (f'--overwrite ' if overwrite else '')
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=[combiner_job],
            job_name='Sample QC hard filters',
        )
    else:
        sample_qc_hardfilter_job = b.new_job(f'Sample QC hard filters [reuse]')

    mt_union_hgdp_for_pca_path = join(
        sample_qc_bucket, 'ancestry', 'mt_union_hgdp_for_pca_path.mt'
    )
    mt_for_pca_path = join(sample_qc_bucket, 'ancestry', 'mt_for_pca.mt')
    provided_pop_ht_path = join(sample_qc_bucket, 'ancestry', 'provided_pop.ht')
    if not can_reuse(
        [
            mt_for_pca_path,
            mt_union_hgdp_for_pca_path,
            provided_pop_ht_path,
        ],
        overwrite,
    ):
        subset_for_pca_job = dataproc.hail_dataproc_job(
            b,
            f'{utils.SCRIPTS_DIR}/sample_qc_subset_mt_for_pca.py '
            + (f'--overwrite ' if overwrite else '')
            + f'--mt {mt_path} '
            f'--input-metadata-ht {input_metadata_ht_path} '
            + (
                f'--pre-computed-hgdp-union-mt {pre_computed_hgdp_unuon_mt_path} '
                if pre_computed_hgdp_unuon_mt_path
                else ''
            )
            + f'--out-hgdp-union-mt {mt_union_hgdp_for_pca_path} '
            f'--out-provided-pop-ht {provided_pop_ht_path} '
            f'--out-mt {mt_for_pca_path} '
            + ('--is-test ' if is_test else '')
            + (f'--pop {pca_pop} ' if pca_pop else '')
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=[combiner_job],
            job_name='Subset MT for PCA',
        )
    else:
        subset_for_pca_job = b.new_job(f'Subset MT for PCA [reuse]')

    relatedness_ht_path = join(out_analysis_bucket, 'relatedness.ht')
    if not can_reuse(relatedness_ht_path, overwrite):
        pcrelate_job = dataproc.hail_dataproc_job(
            b,
            f'{utils.SCRIPTS_DIR}/sample_qc_pcrelate.py '
            + (f'--overwrite ' if overwrite else '')
            + f'--pca-mt {mt_for_pca_path} '
            f'--out-relatedness-ht {relatedness_ht_path} '
            f'--tmp-bucket {sample_qc_bucket}/tmp '
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            # Spark would have problems shuffling on pre-emptible workers
            # (throw SparkException: Job aborted due to stage failure: ShuffleMapStage)
            # so we use num_workers instead of num_secondary_workers here
            num_workers=scatter_count,
            depends_on=[subset_for_pca_job],
            job_name='Run pcrelate',
        )
    else:
        pcrelate_job = b.new_job(f'Run pcrelate [reuse]')

    intermediate_related_samples_to_drop_ht_path = join(
        sample_qc_bucket, 'tmp', 'intermediate_related_samples_to_drop.ht'
    )
    if not can_reuse(intermediate_related_samples_to_drop_ht_path, overwrite):
        cutoffs_d = utils.get_filter_cutoffs(filter_cutoffs_path)
        flag_related_job = dataproc.hail_dataproc_job(
            b,
            f'{utils.SCRIPTS_DIR}/sample_qc_flag_related.py '
            f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
            f'--sex-ht {sex_ht_path} '
            f'--relatedness-ht {relatedness_ht_path} '
            f'--out-ht {intermediate_related_samples_to_drop_ht_path} '
            f'--tmp-bucket {sample_qc_bucket}/tmp '
            f'--max-kin {cutoffs_d["pca"]["max_kin"]} '
            + (f'--overwrite ' if overwrite else '')
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=[sample_qc_hardfilter_job, pcrelate_job],
            job_name='Sample QC flag related',
        )
    else:
        flag_related_job = b.new_job(f'Sample QC flag related [reuse]')

    pca_job, pca_scores_ht_path = _add_ancestry_jobs(
        b=b,
        mt_union_hgdp_path=mt_union_hgdp_for_pca_path,
        provided_pop_ht_path=provided_pop_ht_path,
        num_ancestry_pcs=num_ancestry_pcs,
        tmp_bucket=join(sample_qc_bucket, 'tmp'),
        out_analysis_bucket=out_analysis_bucket,
        out_web_bucket=out_web_bucket,
        overwrite=overwrite,
        scatter_count=scatter_count,
        depends_on=[flag_related_job],
        related_samples_to_drop_ht_path=intermediate_related_samples_to_drop_ht_path,
        billing_project=billing_project,
    )

    inferred_pop_ht_path = join(sample_qc_bucket, 'ancestry', 'inferred_pop.ht')
    regressed_metrics_ht_path = join(sample_qc_bucket, 'regressed_metrics.ht')
    if not can_reuse(
        [
            inferred_pop_ht_path,
            regressed_metrics_ht_path,
        ],
        overwrite,
    ):
        regressed_filters_job = dataproc.hail_dataproc_job(
            b,
            f'{utils.SCRIPTS_DIR}/sample_qc_regressed_filters.py '
            f'--pca-scores-ht {pca_scores_ht_path} '
            f'--provided-pop-ht {provided_pop_ht_path} '
            f'--hail-sample-qc-ht {hail_sample_qc_ht_path} '
            f'{filter_cutoffs_param} '
            f'--n-pcs {num_ancestry_pcs} '
            f'--out-regressed-metrics-ht {regressed_metrics_ht_path} '
            f'--out-inferred-pop-ht {inferred_pop_ht_path} '
            f'--tmp-bucket {sample_qc_bucket}/tmp '
            + (f'--overwrite ' if overwrite else '')
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=[flag_related_job, pca_job, sample_qc_hardfilter_job],
            job_name='Sample QC regressed filters',
        )
    else:
        regressed_filters_job = b.new_job(f'Sample QC regressed filters [reuse]')

    meta_ht_path = join(out_analysis_bucket, 'meta.ht')
    meta_tsv_path = join(out_analysis_bucket, 'meta.tsv')
    if not can_reuse(
        [
            meta_ht_path,
            meta_tsv_path,
        ],
        overwrite,
    ):
        if utils.file_exists(age_csv_path):
            age_csv_param = f'--age-csv {age_csv_path} '
        else:
            age_csv_param = ''
        metadata_qc_job = dataproc.hail_dataproc_job(
            b,
            f'{utils.SCRIPTS_DIR}/sample_qc_write_metadata.py '
            f'--input-metadata-ht {input_metadata_ht_path} '
            f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
            f'--sex-ht {sex_ht_path} '
            f'--custom-qc-ht {custom_qc_ht_path} '
            f'--hail-sample-qc-ht {hail_sample_qc_ht_path} '
            f'--regressed-filtes-ht {regressed_metrics_ht_path} '
            f'--relatedness-ht {relatedness_ht_path} '
            f'--pop-ht {inferred_pop_ht_path} '
            f'{age_csv_param}'
            f'--tmp-bucket {sample_qc_bucket}/tmp '
            f'--out-meta-ht {meta_ht_path} '
            f'--out-meta-tsv {meta_tsv_path} '
            + (f'--overwrite ' if overwrite else '')
            + f'{filter_cutoffs_param} '
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            depends_on=[regressed_filters_job],
            job_name=f'Write sample QC metadata',
        )
    else:
        metadata_qc_job = b.new_job(f'Write sample QC metadata [reuse]')

    if pca_pop:
        _add_ancestry_jobs(
            b=b,
            pop=pca_pop,
            mt_union_hgdp_path=mt_union_hgdp_for_pca_path,
            provided_pop_ht_path=provided_pop_ht_path,
            num_ancestry_pcs=num_ancestry_pcs,
            tmp_bucket=join(sample_qc_bucket, 'tmp'),
            out_analysis_bucket=out_analysis_bucket,
            out_web_bucket=out_web_bucket,
            overwrite=overwrite,
            scatter_count=scatter_count,
            depends_on=[flag_related_job],
            related_samples_to_drop_ht_path=intermediate_related_samples_to_drop_ht_path,
            billing_project=billing_project,
        )

    return metadata_qc_job, hard_filtered_samples_ht_path, meta_ht_path


def _add_ancestry_jobs(
    b: hb.Batch,
    mt_union_hgdp_path: str,
    provided_pop_ht_path: str,
    num_ancestry_pcs: int,
    tmp_bucket: str,
    out_analysis_bucket: str,
    out_web_bucket: str,
    overwrite: bool,
    scatter_count: int,
    pop: Optional[str] = None,
    depends_on: Optional[List[Job]] = None,
    related_samples_to_drop_ht_path: Optional[str] = None,
    billing_project: Optional[str] = None,
) -> Tuple[Job, str]:
    """
    Run PCA and make PCA and loadings plots.

    Returns the PCA job, and the path to PCA scores Table
    """
    pop_tag = pop or 'all'

    ancestry_analysis_bucket = join(out_analysis_bucket, 'ancestry', pop_tag)
    ancestry_web_bucket = join(out_web_bucket, 'ancestry', pop_tag)

    eigenvalues_path = join(ancestry_analysis_bucket, f'eigenvalues_{pop_tag}.txt')
    scores_ht_path = join(ancestry_analysis_bucket, f'scores_{pop_tag}.ht')
    loadings_ht_path = join(ancestry_analysis_bucket, f'loadings_{pop_tag}.ht')
    job_name = f'PCA ({pop or "all"})'
    if not can_reuse([eigenvalues_path, scores_ht_path, loadings_ht_path], overwrite):
        pca_job = dataproc.hail_dataproc_job(
            b,
            f'{utils.SCRIPTS_DIR}/ancestry_pca.py '
            + f'--hgdp-union-mt {mt_union_hgdp_path} '
            + (f'--pop {pop} ' if pop else '')
            + f'--n-pcs {num_ancestry_pcs} '
            f'--out-eigenvalues {eigenvalues_path} '
            f'--out-scores-ht {scores_ht_path} '
            f'--out-loadings-ht {loadings_ht_path} '
            + (
                f'--related-samples-to-drop-ht {related_samples_to_drop_ht_path} '
                if related_samples_to_drop_ht_path
                else ''
            )
            + f'--tmp-bucket {tmp_bucket} '
            + (f'--overwrite ' if overwrite else '')
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=depends_on or [],
            job_name=job_name,
        )
    else:
        pca_job = b.new_job(f'{job_name} [reuse]')

    out_path_ptn = join(ancestry_web_bucket, '{scope}_pc{pci}.{ext}')
    paths = []
    for scope in ['study', 'continental_pop', 'subpop', 'loadings']:
        for ext in ['png', 'html']:
            paths.append(
                out_path_ptn.format(scope=scope, pci=num_ancestry_pcs - 1, ext=ext)
            )
    job_name = f'Plot PCA and loadings ({pop_tag})'
    if not can_reuse(paths, overwrite):
        dataproc.hail_dataproc_job(
            b,
            f'{join(utils.SCRIPTS_DIR, "ancestry_plot.py")} '
            f'--eigenvalues {eigenvalues_path} '
            f'--scores-ht {scores_ht_path} '
            f'--loadings-ht {loadings_ht_path} '
            f'--provided-pop-ht {provided_pop_ht_path} '
            + f'--out-path-pattern {out_path_ptn} '
            + (f'--hail-billing {billing_project} ' if billing_project else ''),
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES + ['selenium'],
            init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
            num_secondary_workers=scatter_count,
            depends_on=[pca_job],
            job_name=job_name,
        )
    else:
        b.new_job(f'{job_name} [reuse]')

    return pca_job, scores_ht_path


def find_inputs(
    output_bucket: str,
    overwrite: bool,
    input_projects: List[str],
    is_test: bool,
    skip_samples: Optional[Collection[str]] = None,
) -> pd.DataFrame:
    """
    Find inputs, make a sample data DataFrame, save to a TSV file
    """
    output_tsv_path = join(output_bucket, 'samples.tsv')

    if utils.can_reuse(output_tsv_path, overwrite):
        logger.info(f'Reading already found DB inputs from {output_tsv_path}')
        samples_df = pd.read_csv(output_tsv_path, sep='\t', na_values='NA').set_index(
            's', drop=False
        )
        return samples_df

    logger.info(
        f'Querying samples from the sample-metadata server '
        f'for the projects: {", ".join(input_projects)}'
    )
    samples_df = sm_utils.find_inputs_from_db(
        input_projects,
        is_test=is_test,
        skip_samples=skip_samples,
    )
    samples_df = sm_utils.add_validation_samples(samples_df)
    samples_df.to_csv(output_tsv_path, index=False, sep='\t', na_rep='NA')
    samples_df = samples_df[pd.notnull(samples_df.s)]
    return samples_df


if __name__ == '__main__':
    main()  # pylint: disable=E1120
