#!/usr/bin/env python3

"""
Hail Batch workflow to perform joint calling, sample QC, and variant QC with VQSR and 
random forest methods on a WGS germline callset.

The workflow is parametrised by the access level, the dataset name, 
batch names and the output version.
"""

import os
import subprocess
import tempfile
import traceback
from os.path import join
from typing import List, Optional, Collection
import logging
import click
import pandas as pd
import hailtop.batch as hb
from analysis_runner import dataproc
from sample_metadata.apis import AnalysisApi
from sample_metadata.models import AnalysisModel, AnalysisStatus, AnalysisType

from joint_calling import utils
from joint_calling import sm_utils
from joint_calling.jobs.variant_qc import add_variant_qc_jobs
from joint_calling.jobs.sample_qc import add_sample_qc_jobs
from joint_calling.jobs import pre_combiner, pedigree

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
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
@click.option(
    '--source-tag',
    'source_tag',
    help='Filter analysis entries by this tag to get GVCF inputs',
)
@click.option(
    '--input-tsv',
    'input_tsv_path',
    help='Instead of pulling samples from the sample-metadata, use this TSV directly',
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
    'ped_fpath',
    help='PED file with family information',
    callback=utils.get_validation_callback(ext='ped', must_exist=True),
)
@click.option(
    '--age-file',
    'age_datas',
    multiple=True,
    help='format: --age-file <tsv-or-csv-path>::<sample-col>::<age-col>. '
    'E.g.: --age-file gs://path_to/age.csv::1::2',
    callback=utils.ColumnInFile.callback,
)
@click.option(
    '--reported-sex-file',
    'reported_sex_datas',
    multiple=True,
    help='format: --reported-sex-file <tsv-or-csv-path>::<sample-col>::<sex-col>. '
    'E.g.: --reported-sex-file gs://path_to/reported_sex.tsv::0::1',
    callback=utils.ColumnInFile.callback,
)
@click.option(
    '--filter-cutoffs-file',
    'filter_cutoffs_path',
    callback=utils.get_validation_callback(ext='yaml', must_exist=True),
    help=f'YAML file with filtering cutoffs. '
    f'Default is the file within the package: {utils.get_filter_cutoffs()}',
)
@click.option(
    '--keep-scratch/--remove-scratch',
    'keep_scratch',
    default=True,
    is_flag=True,
)
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
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option(
    '--check-existence/--no-check-existence',
    'check_existence',
    default=False,
    is_flag=True,
)
@click.option(
    '--add-validation-samples/--no-add-validation-samples',
    'add_validation_samples',
    default=True,
    is_flag=True,
)
@click.option(
    '--assume-gvcfs-are-ready',
    '--validate-smdb',
    'assume_gvcfs_are_ready',
    is_flag=True,
    help='Do not validate the existence of `analysis.output` of SMDB GVCF records',
)
@click.option(
    '--release-related/--release-unrelated',
    'release_related',
    default=False,
    is_flag=True,
    help='Whether to keep or remove related samples from the release',
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,too-many-statements
    output_namespace: str,
    analysis_project: str,
    input_projects: List[str],  # pylint: disable=unused-argument
    source_tag: Optional[str],
    input_tsv_path: Optional[str],  # pylint: disable=unused-argument
    output_version: str,
    output_projects: Optional[List[str]],  # pylint: disable=unused-argument
    use_gnarly_genotyper: bool,
    ped_fpath: Optional[str],
    age_datas: Optional[List[utils.ColumnInFile]],
    reported_sex_datas: Optional[List[utils.ColumnInFile]],
    filter_cutoffs_path: str,
    keep_scratch: bool,
    reuse_scratch_run_id: str,  # pylint: disable=unused-argument
    overwrite: bool,
    skip_samples: Collection[str],
    scatter_count: int,
    pca_pop: Optional[str],
    num_ancestry_pcs: int,
    dry_run: bool,
    check_existence: bool,
    add_validation_samples: bool,
    assume_gvcfs_are_ready: bool,
    release_related: bool,
):  # pylint: disable=missing-function-docstring
    # Determine bucket paths
    if output_namespace in ['test', 'tmp']:
        tmp_bucket_suffix = 'test-tmp'
        project_output_suffix = 'test'  # from crams and gvcfs
        input_sm_projects = [p + '-test' for p in input_projects]
    else:
        tmp_bucket_suffix = 'main-tmp'
        project_output_suffix = 'main'
        input_sm_projects = input_projects

    if output_namespace in ['test', 'main']:
        output_suffix = output_namespace
        output_analysis_suffix = f'{output_namespace}-analysis'
        web_bucket_suffix = f'{output_namespace}-web'
    else:
        output_suffix = 'test-tmp'
        output_analysis_suffix = 'test-tmp/analysis'
        web_bucket_suffix = 'test-tmp/web'

    ptrn = f'gs://cpg-{analysis_project}-{{suffix}}/joint-calling/{output_version}'
    tmp_bucket = ptrn.format(suffix=tmp_bucket_suffix)
    analysis_bucket = ptrn.format(suffix=output_analysis_suffix)
    web_bucket = ptrn.format(suffix=web_bucket_suffix)

    combiner_bucket = f'{tmp_bucket}/combiner'
    raw_combined_mt_path = f'{combiner_bucket}/{output_version}-raw.mt'
    output_bucket = f'gs://cpg-{analysis_project}-{output_suffix}'
    filtered_combined_mt_path = f'{output_bucket}/mt/{output_version}.mt'
    filtered_vcf_ptrn_path = f'{output_bucket}/vcf/{output_version}.chr{{CHROM}}.vcf.gz'
    _create_mt_readme(
        raw_combined_mt_path,
        filtered_combined_mt_path,
        output_bucket,
        filtered_vcf_ptrn_path,
    )

    # Initialize Hail Batch
    billing_project = os.getenv('HAIL_BILLING_PROJECT') or analysis_project
    hail_bucket = os.environ.get('HAIL_BUCKET')
    if not hail_bucket or keep_scratch or reuse_scratch_run_id:
        # Scratch files are large, so we want to use the temporary bucket for them
        hail_bucket = f'{tmp_bucket}/hail'
    logger.info(
        f'Starting hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=hail_bucket.replace('gs://', ''),
        token=os.getenv('HAIL_TOKEN'),
    )
    b = hb.Batch(
        f'Joint calling: {analysis_project}'
        f', version: {output_version}'
        f', data from: {input_tsv_path or ", ".join(input_projects)}'
        + (
            f', limit output projects to {", ".join(output_projects)}'
            if output_projects
            else ''
        )
        + f', namespace: {output_namespace}',
        backend=backend,
    )

    samples_df = find_inputs(
        tmp_bucket=tmp_bucket,
        overwrite=overwrite,
        input_projects=input_sm_projects,
        source_tag=source_tag,
        input_tsv_path=input_tsv_path,
        skip_samples=skip_samples,
        check_existence=check_existence,
        age_datas=age_datas,
        reported_sex_datas=reported_sex_datas,
        add_validation_samples=add_validation_samples,
    )

    (
        samples_df,
        samples_tsv_path,
        pre_combiner_jobs,
    ) = pre_combiner.add_pre_combiner_jobs(
        b=b,
        samples_df=samples_df,
        pre_combiner_bucket=join(tmp_bucket, 'pre_combiner'),
        output_suffix=project_output_suffix,
        overwrite=overwrite,
        analysis_project=analysis_project,
        assume_gvcfs_are_ready=assume_gvcfs_are_ready,
    )

    relatedness_bucket = join(analysis_bucket, 'relatedness')
    somalier_j, ped_fpath, _, somalier_pairs_path = pedigree.pedigree_checks(
        b,
        samples_df=samples_df,
        overwrite=overwrite,
        ped_fpath=ped_fpath,
        output_suffix=project_output_suffix,
        relatedness_bucket=relatedness_bucket,
        web_bucket=join(web_bucket, 'somalier'),
        web_url=f'https://{output_namespace}-web.populationgenomics.org.au/{analysis_project}',
        tmp_bucket=join(tmp_bucket, 'somalier'),
        depends_on=pre_combiner_jobs,
        assume_files_exist=assume_gvcfs_are_ready,
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
        if not utils.can_reuse(raw_combined_mt_path, overwrite):
            combiner_job = dataproc.hail_dataproc_job(
                b,
                f'{utils.SCRIPTS_DIR}/combine_gvcfs.py '
                f'--meta-csv {samples_tsv_path} '
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

    sample_qc_job, hard_filter_ht_path, meta_ht_path = add_sample_qc_jobs(
        b=b,
        mt_path=raw_combined_mt_path,
        samples_tsv_path=samples_tsv_path,
        sample_qc_bucket=join(analysis_bucket, 'sample_qc'),
        ancestry_bucket=join(analysis_bucket, 'ancestry'),
        tmp_bucket=tmp_bucket,
        analysis_bucket=analysis_bucket,
        relatedness_bucket=relatedness_bucket,
        release_related=release_related,
        web_bucket=web_bucket,
        filter_cutoffs_path=filter_cutoffs_path,
        overwrite=overwrite,
        scatter_count=scatter_count,
        combiner_job=combiner_job,
        billing_project=billing_project,
        num_ancestry_pcs=num_ancestry_pcs,
        pca_pop=pca_pop,
        is_test=output_namespace in ['test', 'tmp'],
        somalier_pairs_path=somalier_pairs_path,
        somalier_job=somalier_j,
    )

    var_qc_jobs = add_variant_qc_jobs(
        b=b,
        work_bucket=join(analysis_bucket, 'variant_qc'),
        web_bucket=join(web_bucket, 'variant_qc'),
        raw_combined_mt_path=raw_combined_mt_path,
        hard_filter_ht_path=hard_filter_ht_path,
        meta_ht_path=meta_ht_path,
        out_filtered_combined_mt_path=filtered_combined_mt_path,
        out_filtered_vcf_ptrn_path=filtered_vcf_ptrn_path,
        sample_count=len(samples_df),
        ped_file=ped_fpath,
        overwrite=overwrite,
        vqsr_params_d=utils.get_filter_cutoffs(filter_cutoffs_path)['vqsr'],
        scatter_count=scatter_count,
        is_test=output_namespace in ['test', 'tmp'],
        depends_on=[combiner_job, sample_qc_job],
    )

    # Interacting with the sample metadata server.
    aapi = AnalysisApi()
    # 1. Create a "queued" analysis
    am = AnalysisModel(
        type=AnalysisType('joint-calling'),
        status=AnalysisStatus('queued'),
        output=filtered_combined_mt_path,
        sample_ids=list(samples_df.s),
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
        sm_completed_j.depends_on(*var_qc_jobs)
        logger.info(f'Queueing {am.type} with analysis ID: {aid}')
    except Exception:  # pylint: disable=broad-except
        print(traceback.format_exc())

    b.run(dry_run=dry_run, delete_scratch_on_exit=not keep_scratch, wait=False)


def _create_mt_readme(
    raw_combined_mt_path: str,
    filtered_combined_mt_path: str,
    output_bucket: str,
    filtered_vcf_ptrn_path: str,
):
    local_tmp_dir = tempfile.mkdtemp()
    readme_local_fpath = join(local_tmp_dir, 'README.txt')
    readme_gcs_fpath = join(output_bucket, 'README.txt')
    with open(readme_local_fpath, 'w') as f:
        f.write(f"""\
Unfiltered:
{raw_combined_mt_path}
Soft-filtered matrix table (use `mt.filter_rows(hl.is_missing(mt.filters)` to hard-filter):
{filtered_combined_mt_path}
Hard-filtered per-chromosome VCFs:
{filtered_vcf_ptrn_path}
""")
    subprocess.run(['gsutil', 'cp', readme_local_fpath, readme_gcs_fpath], check=False)


def find_inputs(
    tmp_bucket: str,
    overwrite: bool,
    input_projects: List[str],
    source_tag: Optional[str] = None,
    input_tsv_path: Optional[str] = None,
    skip_samples: Optional[Collection[str]] = None,
    check_existence: bool = True,
    age_datas: Optional[List[utils.ColumnInFile]] = None,
    reported_sex_datas: Optional[List[utils.ColumnInFile]] = None,
    add_validation_samples: bool = True,
) -> pd.DataFrame:
    """
    Find inputs, make a sample data DataFrame, save to a TSV file.

    To specify where to search raw GVCFs, specify `raw_gvcf_source_tag` which
    corresponds with meta.source value in 'in-progress 'analysis entries
    """
    output_tsv_path = input_tsv_path or join(tmp_bucket, 'samples.tsv')

    if utils.can_reuse(output_tsv_path, overwrite):
        logger.info(f'Reading already found DB inputs from {output_tsv_path}')
        samples_df = pd.read_csv(output_tsv_path, sep='\t', na_values='NA').set_index(
            's', drop=False
        )
    else:
        logger.info(
            f'Querying samples from the sample-metadata server '
            f'for the projects: {", ".join(input_projects)}'
        )
        samples_df = sm_utils.find_inputs_from_db(
            input_projects,
            skip_samples=skip_samples,
            check_existence=check_existence,
            source_tag=source_tag,
        )
    if age_datas:
        for age_data in age_datas:
            samples_df = _add_age(samples_df, age_data)
    if reported_sex_datas:
        for reported_sex_data in reported_sex_datas:
            samples_df = _add_reported_sex(samples_df, reported_sex_data)
    if add_validation_samples:
        samples_df = sm_utils.add_validation_samples(samples_df)
    samples_df.to_csv(output_tsv_path, index=False, sep='\t', na_rep='NA')
    samples_df = samples_df[pd.notnull(samples_df.s)]
    return samples_df


def _add_age(samples_df: pd.DataFrame, data: utils.ColumnInFile) -> pd.DataFrame:
    data_by_external_id = data.parse(list(samples_df['external_id']))
    for external_id, value in data_by_external_id.items():
        try:
            age = float(value)
        except ValueError:
            pass
        else:
            samples_df = samples_df.set_index('external_id', drop=False)
            try:
                samples_df.loc[external_id, ['age']] = int(age)
            except ValueError:
                pass
            samples_df = samples_df.set_index('s', drop=False)
    return samples_df


def _add_reported_sex(
    samples_df: pd.DataFrame, data: utils.ColumnInFile
) -> pd.DataFrame:
    data_by_external_id = data.parse(list(samples_df['external_id']))
    for external_id, value in data_by_external_id.items():
        if value in ['M', 'Male', '1']:
            sex = '1'
        elif value in ['F', 'Female', '2']:
            sex = '2'
        else:
            sex = '0'
        samples_df = samples_df.set_index('external_id', drop=False)
        samples_df.loc[external_id, ['sex']] = sex
        samples_df = samples_df.set_index('s', drop=False)
    return samples_df


if __name__ == '__main__':
    main()  # pylint: disable=E1120
