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
from os.path import join, dirname, abspath
from typing import List, Optional, Tuple
import logging
import click
import pandas as pd
import hailtop.batch as hb
from hailtop.batch.job import Job
from analysis_runner import dataproc

from joint_calling import utils
from joint_calling.variant_qc import add_variant_qc_jobs

logger = logging.getLogger('joint-calling')
logger.setLevel('INFO')


@click.command()
@click.option('--callset', '--dataset', 'callset_name', type=str, required=True)
@click.option('--version', 'callset_version', type=str, required=True)
@click.option('--batch', 'callset_batches', type=str, multiple=True)
@click.option(
    '--from',
    'input_namespace',
    type=click.Choice(['main', 'test']),
    help='The bucket namespace to read the inputs from',
)
@click.option(
    '--to',
    'output_namespace',
    type=click.Choice(['main', 'test', 'tmp']),
    help='The bucket namespace to write the results to',
)
@click.option(
    '--existing-mt',
    'existing_mt_path',
    callback=utils.get_validation_callback(ext='mt', must_exist=True),
    help='Path to an existing sparse matrix table to combine with the new data',
)
@click.option(
    '--ped-file',
    'ped_file',
    help='PED file with family information',
    type=str,
    callback=utils.get_validation_callback(ext='ped', must_exist=True),
)
@click.option('--skip-input-meta', 'skip_input_meta', is_flag=True)
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
    '--scatter-count',
    'scatter_count',
    type=int,
    required=True,
    help='Number of secondary workers for Dataproc clusters, as well as the'
    'number of shards to parition data for the AS-VQSR analysis',
)
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option('--run-vqsr/--skip-vqsr', 'run_vqsr', is_flag=True, default=True)
@click.option('--run-rf/--skip-rf', 'run_rf', is_flag=True, default=False)
def main(  # pylint: disable=too-many-arguments,too-many-locals,too-many-statements
    callset_name: str,
    callset_version: str,
    callset_batches: List[str],
    input_namespace: str,
    output_namespace: str,
    existing_mt_path: str,
    ped_file: str,
    skip_input_meta: bool,
    filter_cutoffs_path: str,
    keep_scratch: bool,
    reuse_scratch_run_id: str,  # pylint: disable=unused-argument
    overwrite: bool,
    scatter_count: int,
    dry_run: bool,
    run_vqsr: bool,
    run_rf: bool,
):
    """
    Drive a Hail Batch workflow that creates and submits jobs. A job usually runs
    either a Hail Query script from the scripts folder in this repo using a Dataproc
    cluster; or a GATK command using the GATK or Gnarly image.
    """
    logger.info(f'Enable random forest: {run_rf}')
    logger.info(f'Enable VQSR: {run_vqsr}')

    billing_project = os.getenv('HAIL_BILLING_PROJECT') or callset_name

    if not callset_batches:
        raise click.BadParameter(
            'Please, specify batch numbers with --batch '
            '(can put multiple times, e.g. --batch batch1 --batch batch2)'
        )

    if output_namespace in ['test', 'tmp']:
        tmp_bucket_suffix = 'test-tmp'
    else:
        tmp_bucket_suffix = 'main-tmp'

    work_bucket = (
        f'gs://cpg-{callset_name}-{tmp_bucket_suffix}/joint-calling/{callset_version}'
    )
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
        f'Joint calling: {callset_name}'
        f', version: {callset_version}'
        f', batches: {", ".join(callset_batches)}'
        f', from: {input_namespace}'
        f', to: {output_namespace}',
        backend=backend,
    )
    scripts_dir = abspath(join(dirname(__file__), '..', 'scripts'))

    if output_namespace in ['test', 'main']:
        output_suffix = output_namespace
        output_metadata_suffix = f'{output_namespace}-metadata'
        web_bucket_suffix = f'{output_namespace}-web'
    else:
        output_suffix = 'test-tmp'
        output_metadata_suffix = 'test-tmp'
        web_bucket_suffix = 'test-tmp'

    output_metadata_bucket = f'gs://cpg-{callset_name}-{output_metadata_suffix}/joint-calling/{callset_version}'
    web_bucket = (
        f'gs://cpg-{callset_name}-{web_bucket_suffix}/joint-calling/{callset_version}'
    )
    mt_output_bucket = f'gs://cpg-{callset_name}-{output_suffix}/mt'

    raw_combined_mt_path = f'{mt_output_bucket}/{callset_version}-raw.mt'
    # pylint: disable=unused-variable
    filtered_combined_mt_path = f'{mt_output_bucket}/{callset_version}.mt'

    combiner_bucket = f'{work_bucket}/combiner'
    sample_qc_bucket = f'{work_bucket}/sample_qc'

    filter_cutoffs_d = utils.get_filter_cutoffs(filter_cutoffs_path)

    samples_df, samples_csv_path, pre_combiner_jobs = _add_pre_combiner_jobs(
        b=b,
        input_gvcfs_bucket=f'gs://cpg-{callset_name}-{input_namespace}/gvcf',
        input_metadata_bucket=f'gs://cpg-{callset_name}-{input_namespace}-metadata'
        if not skip_input_meta
        else None,
        work_bucket=join(work_bucket, 'pre-combine'),
        output_bucket=combiner_bucket,
        callset_batches=callset_batches,
        overwrite=overwrite,
    )

    if overwrite or not utils.file_exists(raw_combined_mt_path):
        combiner_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/combine_gvcfs.py '
            f'--meta-csv {samples_csv_path} '
            f'--out-mt {raw_combined_mt_path} '
            f'--bucket {combiner_bucket}/work '
            + (f'--existing-mt {existing_mt_path} ' if existing_mt_path else '')
            + f'--hail-billing {billing_project} '
            f'--n-partitions {scatter_count * 25}',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=pre_combiner_jobs,
            job_name='Combine GVCFs',
        )
    else:
        combiner_job = b.new_job('Combine GVCFs [reuse]')

    info_ht_path = join(sample_qc_bucket, 'info.ht')
    info_split_ht_path = join(sample_qc_bucket, 'info-split.ht')
    info_vcf_path = join(sample_qc_bucket, 'info.vcf')
    if overwrite or any(
        not utils.file_exists(fp)
        for fp in [info_ht_path, info_split_ht_path, info_vcf_path]
    ):
        generate_info_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/generate_info_ht.py --overwrite '
            f'--mt {raw_combined_mt_path} '
            f'--out-info-ht {info_ht_path} '
            f'--out-split-info-ht {info_split_ht_path} '
            f'--out-info-vcf {info_vcf_path}',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            # Adding more workers as this is a much longer step
            num_secondary_workers=scatter_count,
            depends_on=[combiner_job],
            job_name='Generate info',
        )
    else:
        generate_info_job = b.new_job('Generate info [reuse]')

    hard_filtered_samples_ht_path = join(sample_qc_bucket, 'hard_filters.ht')
    meta_ht_path = join(output_metadata_bucket, 'meta.ht')
    meta_tsv_path = join(output_metadata_bucket, 'meta.tsv')
    relatedness_ht_path = join(output_metadata_bucket, 'relatedness.ht')
    if overwrite or any(
        not utils.file_exists(fp)
        for fp in [hard_filtered_samples_ht_path, meta_ht_path]
    ):
        age_csv = f'gs://cpg-{callset_name}-{input_namespace}-metadata/age.csv'
        if utils.file_exists(age_csv):
            age_csv_param = f'--age-csv {age_csv} '
        else:
            age_csv_param = ''

        if filter_cutoffs_path:
            gcs_path = join(work_bucket, 'filter-cutoffs.yaml')
            subprocess.run(['gsutil', 'cp', filter_cutoffs_path, gcs_path], check=False)
            filter_cutoffs_param = f'--filter-cutoffs-file {gcs_path}'
        else:
            filter_cutoffs_param = ''

        sample_qc_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/sample_qc.py {filter_cutoffs_param} --overwrite '
            f'--mt {raw_combined_mt_path} '
            f'--info-ht {info_ht_path} '
            f'{age_csv_param}'
            f'--meta-csv {samples_csv_path} '
            f'--bucket {combiner_bucket} '
            f'--out-hardfiltered-samples-ht {hard_filtered_samples_ht_path} '
            f'--out-meta-ht {meta_ht_path} '
            f'--out-meta-tsv {meta_tsv_path} '
            f'--out-relatedness-ht {relatedness_ht_path} '
            f'--hail-billing {billing_project} ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=[combiner_job, generate_info_job],
            job_name='Sample QC',
        )
    else:
        sample_qc_job = b.new_job('Sample QC [reuse]')

    if run_rf or run_vqsr:
        var_qc_job, var_qc_final_filter_ht = add_variant_qc_jobs(
            b=b,
            work_bucket=join(work_bucket, 'variant_qc'),
            analysis_bucket=join(web_bucket, 'variant_qc'),
            raw_combined_mt_path=raw_combined_mt_path,
            info_split_ht_path=info_split_ht_path,
            hard_filter_ht_path=hard_filtered_samples_ht_path,
            meta_ht_path=meta_ht_path,
            samples_df=samples_df,
            sample_qc_job=sample_qc_job,
            scripts_dir=scripts_dir,
            ped_file=ped_file,
            overwrite=overwrite,
            run_rf=run_rf,
            vqsr_params_d=filter_cutoffs_d['vqsr'],
            scatter_count=scatter_count,
        )
        if overwrite or not utils.file_exists(filtered_combined_mt_path):
            finalised_mt_job = dataproc.hail_dataproc_job(
                b,
                f'{scripts_dir}/make_finalised_mt.py --overwrite '
                f'--mt {raw_combined_mt_path} '
                f'--var-qc-final-filter-ht {var_qc_final_filter_ht} '
                f'--out-mt {filtered_combined_mt_path} '
                f'--meta-ht {meta_ht_path} ',
                max_age='8h',
                packages=utils.DATAPROC_PACKAGES,
                num_secondary_workers=scatter_count,
                depends_on=[var_qc_job],
                job_name='Making final MT',
            )
        else:
            finalised_mt_job = b.new_job('Making final MT [reuse]')

    else:
        var_qc_job = b.new_job('Var QC [skip]')

    b.run(dry_run=dry_run, delete_scratch_on_exit=not keep_scratch)


def _add_pre_combiner_jobs(
    b: hb.Batch,
    input_gvcfs_bucket: str,
    input_metadata_bucket: Optional[str],
    work_bucket: str,
    output_bucket: str,
    callset_batches: List[str],
    overwrite: bool,
) -> Tuple[pd.DataFrame, str, List[Job]]:
    """
    Add jobs that prepare GVCFs for the combiner, if needed.
    :param input_gvcfs_bucket: bucket with GVCFs batches as subfolders
    :param input_metadata_bucket: bucket with metadata batches as subfolders
    :param work_bucket: bucket to write intermediate files to
    :param output_bucket: bucket to write the GVCF combiner inputs to
    :param callset_batches: list of the dataset batches identifiers
    (where "batch" refers to just a set of sampels, and not to be confused with
    Hail Batch service)
    :param overwrite: ignore existing intermediate files
    :return: a Tuple of: a pandas dataframe with the sample metadata, a CSV file
    corresponding to that dataframe, and a list of jobs to wait for before
    submitting the combiner job
    """

    # File with the pointers to GVCFs to process along with metdata.
    # If it doesn't exist, we trigger a utuils.find_inputs(combiner_bucket) function
    # to find the GVCFs and the metadata given the requested batch ids.
    input_samples_csv_path = join(work_bucket, 'samples.csv')
    # Raw GVCFs need pre-processing before passing to the combiner. If the following
    # file exists, we assume the samples are pre-processed; otherwise, we add Batch
    # jobs to do the pre-processing.
    combiner_ready_samples_csv_path = join(output_bucket, 'samples.csv')
    subset_gvcf_jobs = []
    if not overwrite and utils.file_exists(combiner_ready_samples_csv_path):
        samples_df = pd.read_csv(combiner_ready_samples_csv_path, sep='\t').set_index(
            's', drop=False
        )
        samples_df = samples_df[pd.notnull(samples_df.s)]
    else:
        if not overwrite and utils.file_exists(input_samples_csv_path):
            samples_df = pd.read_csv(input_samples_csv_path, sep='\t').set_index(
                's', drop=False
            )
        else:
            input_gvcf_buckets = []
            input_metadata_buckets = []
            for cb in callset_batches:
                input_gvcf_buckets.append(join(input_gvcfs_bucket, cb))
                if input_metadata_bucket:
                    input_metadata_buckets.append(join(input_metadata_bucket, cb))
            samples_df = utils.find_inputs(input_gvcf_buckets, input_metadata_buckets)

        samples_df = samples_df[pd.notnull(samples_df.s)]
        gvcfs = [
            b.read_input_group(**{'g.vcf.gz': gvcf, 'g.vcf.gz.tbi': gvcf + '.tbi'})
            for gvcf in list(samples_df.gvcf)
        ]
        subset_gvcf_jobs = _add_prep_gvcfs_for_combiner_steps(
            b=b,
            gvcfs=gvcfs,
            samples_df=samples_df,
            output_gvcf_bucket=join(output_bucket, 'gvcf'),
            overwrite=overwrite,
        )
        samples_df.to_csv(
            combiner_ready_samples_csv_path, index=False, sep='\t', na_rep='NA'
        )
        logger.info(
            f'Saved metadata with updated GVCFs to '
            f'{combiner_ready_samples_csv_path}'
        )

    return samples_df, combiner_ready_samples_csv_path, subset_gvcf_jobs


def _add_prep_gvcfs_for_combiner_steps(
    b,
    gvcfs,
    samples_df: pd.DataFrame,
    output_gvcf_bucket: str,
    overwrite: bool,
    reuse_scratch_run_id: Optional[str] = None,
    hail_bucket: Optional[str] = None,
) -> List[Job]:
    """
    Add steps required to prepare GVCFs from combining
    :param b:
    :param gvcfs:
    :param samples_df: pandas dataframe with metadata and raw GVCF paths
    :param output_gvcf_bucket: bucket to write the combiner-ready GVCFs to
    :param overwrite: overwrite existing intemridate files
    :param reuse_scratch_run_id: ID of a Batch run to reuse the intermediate files
    :param hail_bucket: bucket to find the previous Batch run intermediate files
    :return: list of Batch jobs
    """

    # If we want to re-use the intermediate files from a previous run,
    # find them using the Batch run ID of a previous run with --keep-scratch enabled.
    found_reblocked_gvcf_paths = [
        join(
            str(hail_bucket),
            'batch',
            reuse_scratch_run_id,
            str(job_num),
            'output_gvcf.g.vcf.gz',
        )
        if reuse_scratch_run_id
        else None
        for job_num in range(1, 1 + len(gvcfs))
    ]
    reblocked_gvcfs = [
        b.read_input_group(
            **{
                'vcf.gz': found_gvcf_path,
                'vcf.gz.tbi': found_gvcf_path + '.tbi',
            }
        )
        if (not overwrite and found_gvcf_path and utils.file_exists(found_gvcf_path))
        else _add_reblock_gvcfs_step(b, input_gvcf).output_gvcf
        for found_gvcf_path, input_gvcf in zip(found_reblocked_gvcf_paths, gvcfs)
    ]

    noalt_regions = b.read_input('gs://cpg-reference/hg38/v0/noalt.bed')
    subset_gvcf_jobs = [
        _add_subset_noalt_step(
            b,
            input_gvcf=input_gvcf,
            output_gvcf_path=output_gvcf_path,
            noalt_regions=noalt_regions,
        )
        if not utils.file_exists(output_gvcf_path) or overwrite
        else b.new_job('SubsetToNoalt [reuse]')
        for input_gvcf, output_gvcf_path in [
            (gvcf, join(output_gvcf_bucket, f'{sample}.g.vcf.gz'))
            for sample, gvcf in zip(list(samples_df.s), reblocked_gvcfs)
        ]
    ]
    for sn in samples_df.s:
        samples_df.loc[sn, ['gvcf']] = join(output_gvcf_bucket, sn + '.g.vcf.gz')
    return subset_gvcf_jobs


def _add_reblock_gvcfs_step(
    b: hb.Batch,
    input_gvcf: hb.ResourceGroup,
) -> Job:
    """
    Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration
    """
    j = b.new_job('ReblockGVCF')
    j.image(utils.GATK_DOCKER)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'30G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )

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
    noalt_regions: str,
) -> Job:
    """
    1. Subset GVCF to main chromosomes to avoid downstream errors
    2. Removes the DS INFO field that is added to some HGDP GVCFs to avoid errors
       from Hail about mismatched INFO annotations
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
    j.command(
        f"""set -e

    bcftools view \\
        {input_gvcf['g.vcf.gz']} \\
        -T {noalt_regions} \\
        | bcftools annotate -x INFO/DS \\
        -o {j.output_gvcf['g.vcf.gz']} \\
        -Oz

    bcftools index --tbi {j.output_gvcf['g.vcf.gz']}
        """
    )
    b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
    return j


if __name__ == '__main__':
    main()  # pylint: disable=E1120
