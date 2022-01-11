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
from os.path import join
from typing import List, Optional, Collection
import logging
import click
import pandas as pd
import hailtop.batch as hb

from joint_calling.dataproc import add_job
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
    default=False,
    help='if an intermediate or a final file exists, skip running the code '
    'that generates it.',
)
@click.option(
    '--skip-sample',
    'skip_samples',
    multiple=True,
    help='Don\'t process specified samples. Can be set multiple times.',
)
@click.option(
    '--force-sample',
    'force_samples',
    multiple=True,
)
@click.option(
    '--scatter-count',
    'scatter_count',
    type=int,
    required=True,
    help='Number of secondary workers for Dataproc clusters, as well as the'
    'number of shards to parition data for the AS-VQSR analysis',
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
    '--do-not-query-smdb-for-gvcfs',
    'do_not_query_smdb_for_gvcfs',
    is_flag=True,
)
@click.option(
    '--assume-gvcfs-exist',
    'assume_gvcfs_exist',
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
@click.option(
    '--skip-somalier',
    'skip_somalier',
    default=False,
    is_flag=True,
    help='Do not generate somaleir fingerprints. Use pc_relate',
)
@click.option(
    '--combiner-branch-factor',
    'combiner_branch_factor',
    type=click.INT,
    help='Branch factor for VCF combiner (see https://hail.is/docs/0.2/experimental/vcf_combiner.html#pain-points)'
)
@click.option(
    '--combiner-batch-size',
    'combiner_batch_size',
    type=click.INT,
    help='Size of a batch for phase1 of VCF combiner (see https://hail.is/docs/0.2/experimental/vcf_combiner.html#pain-points)'
)
@click.option(
    '--highmem-workers',
    'highmem_workers',
    is_flag=True,
    help='Use highmem workers in dataproc'
)
@click.option(
    '--skip-missing-qc',
    'skip_missing_qc',
    is_flag=True,
    help='Drop samples without QC provided'
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,too-many-statements
    output_namespace: str,
    analysis_project: str,
    input_projects: List[str],  # pylint: disable=unused-argument
    source_tag: Optional[str],
    input_tsv_path: Optional[str],  # pylint: disable=unused-argument
    output_version: str,
    output_projects: Optional[List[str]],  # pylint: disable=unused-argument
    ped_fpath: Optional[str],
    age_datas: Optional[List[utils.ColumnInFile]],
    reported_sex_datas: Optional[List[utils.ColumnInFile]],
    filter_cutoffs_path: str,
    keep_scratch: bool,
    reuse_scratch_run_id: str,  # pylint: disable=unused-argument
    overwrite: bool,
    skip_samples: Optional[Collection[str]],
    force_samples: Optional[Collection[str]],
    scatter_count: int,
    dry_run: bool,
    check_existence: bool,
    add_validation_samples: bool,
    do_not_query_smdb_for_gvcfs: bool,
    assume_gvcfs_exist: bool,
    release_related: bool,
    skip_somalier: bool,
    combiner_branch_factor: Optional[int],
    combiner_batch_size: Optional[int],
    highmem_workers: bool,
    skip_missing_qc: bool,
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

    combiner_bucket = f'{analysis_bucket}/combiner'
    raw_combined_mt_path = f'{combiner_bucket}/{output_version}-raw.mt'
    output_bucket = f'gs://cpg-{analysis_project}-{output_suffix}'
    filtered_combined_mt_path = f'{output_bucket}/mt/{output_version}.mt'
    filtered_vcf_ptrn_path = f'{output_bucket}/vcf/{output_version}.chr{{CHROM}}.vcf.bgz'
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
        analysis_bucket=analysis_bucket,
        overwrite=False,
        input_projects=input_sm_projects,
        source_tag=source_tag,
        input_tsv_path=input_tsv_path,
        skip_samples=skip_samples,
        check_existence=check_existence,
        age_datas=age_datas,
        reported_sex_datas=reported_sex_datas,
        add_validation_samples=add_validation_samples,
        do_not_query_smdb_for_gvcfs=do_not_query_smdb_for_gvcfs,
    )

    (
        samples_df,
        samples_tsv_path,
        pre_combiner_jobs,
    ) = pre_combiner.add_pre_combiner_jobs(
        b=b,
        samples_df=samples_df,
        tmp_bucket=tmp_bucket,
        analysis_bucket=analysis_bucket,
        output_suffix=project_output_suffix,
        overwrite=overwrite,
        analysis_project=analysis_project,
        force_samples=force_samples or [],
        skip_samples=skip_samples or [],
        assume_gvcfs_exist=assume_gvcfs_exist,
    )

    relatedness_bucket = join(analysis_bucket, 'relatedness')
    somalier_j = None
    somalier_pairs_path = None
    if not skip_somalier:
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
            assume_files_exist=assume_gvcfs_exist,
        )

    # Combiner takes advantage of autoscaling cluster policies
    # to reduce costs for the work that uses only the driver machine:
    # https://hail.is/docs/0.2/experimental/vcf_combiner.html#pain-points
    if scatter_count > 100:
        autoscaling_workers = '200'
    elif scatter_count > 50:
        autoscaling_workers = '100'
    else:
        autoscaling_workers = '50'

    if not utils.can_reuse(raw_combined_mt_path, overwrite):
        combiner_job = add_job(
            b,
            f'{utils.SCRIPTS_DIR}/combine_gvcfs.py '
            f'--meta-tsv {samples_tsv_path} '
            f'--out-mt {raw_combined_mt_path} '
            f'--bucket {combiner_bucket}/work '
            f'--hail-billing {billing_project} ' +
            (f'--branch-factor {combiner_branch_factor} ' if combiner_branch_factor else '') +
            (f'--batch-size {combiner_batch_size} ' if combiner_batch_size else '') +
            f'--n-partitions {scatter_count * 25}',
            job_name='Combine GVCFs',
            num_workers=0,
            highmem=highmem_workers,
            autoscaling_policy=f'vcf-combiner-{autoscaling_workers}',
            long=True,
        )
    else:
        combiner_job = b.new_job('Combine GVCFs [reuse]')
    combiner_job.depends_on(*pre_combiner_jobs)

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
        is_test=output_namespace in ['test', 'tmp'],
        somalier_pairs_path=somalier_pairs_path,
        somalier_job=somalier_j,
        highmem_workers=highmem_workers,
    )

    add_variant_qc_jobs(
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
        project_name=analysis_project,
        highmem_workers=highmem_workers,
    )

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
    analysis_bucket: str,
    overwrite: bool,
    input_projects: List[str],
    source_tag: Optional[str] = None,
    input_tsv_path: Optional[str] = None,
    skip_samples: Optional[Collection[str]] = None,
    check_existence: bool = True,
    age_datas: Optional[List[utils.ColumnInFile]] = None,
    reported_sex_datas: Optional[List[utils.ColumnInFile]] = None,
    add_validation_samples: bool = True,
    do_not_query_smdb_for_gvcfs: bool = False,
) -> pd.DataFrame:
    """
    Find inputs, make a sample data DataFrame, save to a TSV file.

    To specify where to search raw GVCFs, specify `raw_gvcf_source_tag` which
    corresponds with meta.source value in 'in-progress 'analysis entries
    """
    output_tsv_path = input_tsv_path or join(analysis_bucket, 'samples.tsv')

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
            do_not_query_smdb_for_gvcfs=do_not_query_smdb_for_gvcfs,
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
