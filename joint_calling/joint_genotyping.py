#!/usr/bin/env python3
# pylint: skip-file
# type: ignore
# flake8: noqa

"""
Functions and batch jobs for joint-genotyping
"""

import json
import logging
import subprocess
from os.path import join, dirname, abspath, splitext, basename
from typing import Optional, List, Tuple, Set, Dict
import pandas as pd
import hailtop.batch as hb
from hailtop.batch.job import Job
import utils
import sm_utils
from joint_calling import resources

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def make_joint_genotype_jobs(
    b: hb.Batch,
    genomicsdb_bucket: str,
    samples_df: pd.DataFrame,
    reference: hb.ResourceGroup,
    dbsnp: str,
    out_bucket: str,
    tmp_bucket: str,
    local_tmp_dir: str,
    overwrite: bool,
    depends_on: Optional[List[Job]] = None,
    analysis_project: str = None,
    completed_analysis: Optional[sm_utils.Analysis] = None,
) -> Tuple[Job, str]:
    """
    Assumes all samples have a 'file' of 'type'='gvcf' in `samples_df`.
    Adds samples to the GenomicsDB and runs joint genotyping on them.
    Outputs a multi-sample VCF under `output_vcf_path`.
    """
    is_small_callset = len(samples_df) < 1000
    # 1. For small callsets, we don't apply the ExcessHet filtering.
    # 2. For small callsets, we gather the VCF shards and collect QC metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics
    # collected from them.
    is_huge_callset = len(samples_df) >= 100000
    # For huge callsets, we allocate more memory for the SNPs Create Model step

    genomics_gcs_path_per_interval = dict()
    for idx in range(utils.NUMBER_OF_GENOMICS_DB_INTERVALS):
        genomics_gcs_path_per_interval[idx] = join(
            genomicsdb_bucket,
            f'interval_{idx}_outof_{utils.NUMBER_OF_GENOMICS_DB_INTERVALS}.tar',
        )
    # Determining which samples to add. Using the first interval, so the assumption
    # is that all DBs have the same set of samples.
    (
        sample_names_to_add,
        sample_names_will_be_in_db,
        sample_names_already_added,
        updating_existing_db,
        sample_map_bucket_path,
    ) = _samples_to_add_to_db(
        genomicsdb_gcs_path=genomics_gcs_path_per_interval[0],
        interval_idx=0,
        local_tmp_dir=local_tmp_dir,
        samples=samples,
        tmp_bucket=tmp_bucket,
        gvcf_by_sid=gvcf_by_sid,
    )
    if not sample_names_to_add:
        logger.info('')

    samples_hash = utils.hash_sample_names(sample_names_will_be_in_db)
    gathered_vcf_path = join(
        out_bucket, 'jointly-called', 'tmp', 'gathered', samples_hash + '.vcf.gz'
    )
    job_name = 'Joint-calling+VQSR'
    vqsred_vcf_path = join(out_bucket, 'jointly-called', samples_hash + '.vcf.gz')
    if completed_analysis:
        if utils.file_exists(vqsred_vcf_path):
            logger.info(
                f'Completed analysis and the output for these samples exist: '
                f'"{vqsred_vcf_path}", reusing.'
            )
            j = b.new_job('Joint calling+VQSR [reuse]')
            return j, vqsred_vcf_path
        else:
            logger.warning(
                f'Joint-calling completed analysis for these samples exists, '
                f'but the output file is missing. Setting status to "failed" and '
                f'rerunning the joint-calling analysis.'
            )
            aapi.update_analysis_status(
                completed_analysis.id,
                AnalysisUpdateModel(status='failed'),
            )
    else:
        logger.info(
            f'Completed joint-calling analysis does not exist for this set of samples. '
            f'Submitting the joint-calling and VQSR jobs.'
        )

    am = AnalysisModel(
        sample_ids=[s['id'] for s in samples],
        type='joint-calling',
        output=vqsred_vcf_path,
        status='queued',
    )
    if utils.can_reuse(vqsred_vcf_path, overwrite=overwrite):
        # Create a "completed" analysis and return an empty job
        am.status = 'completed'
        if analysis_project:
            aapi.create_new_analysis(project=analysis_project, analysis_model=am)
        return b.new_job(f'{job_name} [reuse]'), vqsred_vcf_path

    intervals_j = _add_split_intervals_job(
        b=b,
        interval_list=resources.UNPADDED_INTERVALS,
        scatter_count=utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
        ref_fasta=resources.REF_FASTA,
    )

    if analysis_project:
        # Interacting with the sample metadata server:
        # 1. Create a "queued" analysis
        aid = aapi.create_new_analysis(project=analysis_project, analysis_model=am)
        # 2. Queue a job that updates the status to "in-progress"
        sm_in_progress_j = utils.make_sm_in_progress_job(
            b,
            project=analysis_project,
            analysis_id=aid,
            analysis_type='joint-calling',
        )
        # 2. Queue a job that updates the status to "completed"
        sm_completed_j = utils.make_sm_completed_job(
            b,
            project=analysis_project,
            analysis_id=aid,
            analysis_type='joint-calling',
        )
        # Set up dependencies
        intervals_j.depends_on(sm_in_progress_j)
        if depends_on:
            sm_in_progress_j.depends_on(*depends_on)
        logger.info(f'Queueing {am.type} with analysis ID: {aid}')
    else:
        if depends_on:
            intervals_j.depends_on(*depends_on)
        sm_completed_j = None

    import_gvcfs_job_per_interval = dict()
    if sample_names_to_add:
        for idx in range(utils.NUMBER_OF_GENOMICS_DB_INTERVALS):
            import_gvcfs_job, _ = _add_import_gvcfs_job(
                b=b,
                genomicsdb_gcs_path=genomics_gcs_path_per_interval[idx],
                sample_names_to_add=sample_names_to_add,
                sample_names_to_skip=sample_names_already_added,
                sample_names_will_be_in_db=sample_names_will_be_in_db,
                updating_existing_db=updating_existing_db,
                sample_map_bucket_path=sample_map_bucket_path,
                interval=intervals_j.intervals[f'interval_{idx}'],
                interval_idx=idx,
                number_of_intervals=utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
                depends_on=[intervals_j],
            )
            import_gvcfs_job_per_interval[idx] = import_gvcfs_job

    make_site_only_jobs = []
    scattered_vcf_by_interval: Dict[int, hb.ResourceGroup] = dict()

    for idx in range(utils.NUMBER_OF_GENOMICS_DB_INTERVALS):
        samples_hash = utils.hash_sample_names(sample_names_will_be_in_db)
        site_only_vcf_path = join(
            tmp_bucket, 'jointly-called', samples_hash, f'interval_{idx}.vcf.gz'
        )
        if utils.can_reuse(site_only_vcf_path, overwrite):
            make_site_only_jobs.append(b.new_job('Joint genotyping [reuse]'))
            scattered_vcf_by_interval[idx] = b.read_input_group(
                **{
                    'vcf.gz': site_only_vcf_path,
                    'vcf.gz.tbi': site_only_vcf_path + '.tbi',
                }
            )
        else:
            genotype_vcf_job = _add_gnarly_genotyper_job(
                b,
                genomicsdb_path=genomics_gcs_path_per_interval[idx],
                interval=intervals_j.intervals[f'interval_{idx}'],
                reference=reference,
                dbsnp=dbsnp,
                overwrite=overwrite,
                interval_idx=idx,
                number_of_samples=len(sample_names_will_be_in_db),
                number_of_intervals=utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
            )
            if import_gvcfs_job_per_interval.get(idx):
                genotype_vcf_job.depends_on(import_gvcfs_job_per_interval.get(idx))

            vcf = genotype_vcf_job.output_vcf
            if not is_small_callset:
                exccess_filter_job = _add_exccess_het_filter(
                    b,
                    input_vcf=vcf,
                    overwrite=overwrite,
                    interval=intervals_j.intervals[f'interval_{idx}'],
                )
                vcf = exccess_filter_job.output_vcf
            make_site_only_job = _add_make_sites_only_job(
                b,
                input_vcf=vcf,
                overwrite=overwrite,
            )
            make_site_only_jobs.append(make_site_only_job)
            scattered_vcf_by_interval[idx] = make_site_only_job.output_vcf

    final_gathered_vcf_job = _add_final_gather_vcf_job(
        b,
        input_vcfs=list(scattered_vcf_by_interval.values()),
        overwrite=overwrite,
        output_vcf_path=gathered_vcf_path,
    )
    tmp_vqsr_bucket = join(tmp_bucket, 'vqsr', samples_hash)
    vqsr_job = make_vqsr_jobs(
        b,
        input_vcf_gathered=gathered_vcf_path,
        input_vcfs_scattered=list(scattered_vcf_by_interval.values()),
        is_small_callset=is_small_callset,
        is_huge_callset=is_huge_callset,
        work_bucket=tmp_vqsr_bucket,
        web_bucket=tmp_vqsr_bucket,
        depends_on=[final_gathered_vcf_job],
        intervals=intervals_j.intervals,
        scatter_count=utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
        output_vcf_path=vqsred_vcf_path,
    )
    if sm_completed_j:
        sm_completed_j.depends_on(vqsr_job)
        return sm_completed_j, vqsred_vcf_path
    else:
        return vqsr_job, vqsred_vcf_path


def _samples_to_add_to_db(
    genomicsdb_gcs_path,
    interval_idx,
    local_tmp_dir,
    samples,
    tmp_bucket: str,
    gvcf_by_sid: Dict[str, str],
) -> Tuple[Set[str], Set[str], Set[str], bool, str]:
    if utils.file_exists(join(genomicsdb_gcs_path, 'callset.json')):
        # Check if sample exists in the DB already
        genomicsdb_metadata = join(local_tmp_dir, f'callset-{interval_idx}.json')
        # This command will download the DB metadata file locally.
        # The `-O` argument to `tar` means "write the file being extracted to the stdout",
        # and the file to be extracted is specified as a positional argument to `tar`.
        cmd = (
            f'gsutil cat {genomicsdb_gcs_path} | '
            f'tar -O --extract workspace/callset.json > {genomicsdb_metadata}'
        )
        logger.info(cmd)
        subprocess.run(cmd, check=False, shell=True)

        with open(genomicsdb_metadata) as f:
            db_metadata = json.load(f)
        sample_names_in_db = set(s['sample_name'] for s in db_metadata['callsets'])
        sample_names_requested = set([s['id'] for s in samples])
        sample_names_to_add = sample_names_requested - sample_names_in_db
        sample_names_to_remove = sample_names_in_db - sample_names_requested
        if sample_names_to_remove:
            # GenomicsDB doesn't support removing, so creating a new DB
            updating_existing_db = False
            sample_names_already_added = set()
            sample_names_to_add = {s['id'] for s in samples}
            sample_names_will_be_in_db = sample_names_to_add
            logger.info(
                f'GenomicDB {genomicsdb_gcs_path} exists, but '
                f'{len(sample_names_to_remove)} samples need '
                f'to be removed: {", ".join(sample_names_to_remove)}, so creating a new '
                f'DB with {len(sample_names_will_be_in_db)} samples: '
                f'{", ".join(sample_names_will_be_in_db)}'
            )
        else:
            updating_existing_db = True
            sample_names_will_be_in_db = sample_names_in_db | sample_names_to_add
            sample_names_already_added = sample_names_requested & sample_names_in_db
            if sample_names_already_added:
                logger.info(
                    f'{len(sample_names_already_added)} samples '
                    f'{", ".join(sample_names_already_added)} already exist in the DB '
                    f'{genomicsdb_gcs_path}, skipping adding them.'
                )
            if sample_names_to_remove:
                logger.info(
                    f'There are {len(sample_names_to_remove)} samples that need to be '
                    f'removed from the DB {genomicsdb_gcs_path}: '
                    f'{", ".join(sample_names_to_remove)}. Re-creating the DB '
                    f'using the updated set of samples'
                )
            elif sample_names_to_add:
                logger.info(
                    f'Will add {len(sample_names_to_add)} samples '
                    f'{", ".join(sample_names_to_add)} into the DB {genomicsdb_gcs_path}'
                )
            else:
                logger.warning(
                    f'Nothing will be added into the DB {genomicsdb_gcs_path}'
                )
    else:
        # Initiate new DB
        sample_names_already_added = set()
        sample_names_to_add = {s['id'] for s in samples}
        sample_names_will_be_in_db = sample_names_to_add
        updating_existing_db = False
        logger.info(
            f'GenomicDB {genomicsdb_gcs_path} doesn\'t exist, so creating a new one '
            f'with {len(sample_names_to_add)} samples: {", ".join(sample_names_to_add)}'
        )

    sample_map_bucket_path = join(tmp_bucket, 'work', 'sample_name.csv')
    local_sample_map_fpath = join(local_tmp_dir, 'sample_name.csv')
    with open(local_sample_map_fpath, 'w') as f:
        for sid in sample_names_to_add:
            f.write('\t'.join([sid, gvcf_by_sid[sid]]) + '\n')
    subprocess.run(
        f'gsutil cp {local_sample_map_fpath} {sample_map_bucket_path}',
        check=False,
        shell=True,
    )

    return (
        sample_names_to_add,
        sample_names_will_be_in_db,
        sample_names_already_added,
        updating_existing_db,
        sample_map_bucket_path,
    )


def _add_import_gvcfs_job(
    b: hb.Batch,
    genomicsdb_gcs_path: str,
    sample_names_to_add: Set[str],
    sample_names_to_skip: Set[str],
    sample_names_will_be_in_db: Set[str],
    updating_existing_db: bool,
    sample_map_bucket_path: str,
    interval: hb.ResourceFile,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    depends_on: Optional[List[Job]] = None,
) -> Tuple[Optional[Job], Set[str]]:
    """
    Add GVCFs to a genomics database (or create a new instance if it doesn't exist)
    Returns a Job, or None if no new samples to add
    """
    if updating_existing_db:
        # Update existing DB
        genomicsdb_param = f'--genomicsdb-update-workspace-path {genomicsdb_gcs_path}'
        job_name = 'Adding to GenomicsDB'
    else:
        # Initiate new DB
        genomicsdb_param = f'--genomicsdb-workspace-path {genomicsdb_gcs_path}'
        job_name = 'Creating GenomicsDB'

    sample_map = b.read_input(sample_map_bucket_path)

    if interval_idx is not None:
        job_name += f' {interval_idx}/{number_of_intervals}'

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.cpu(16)
    java_mem = 16
    j.memory('lowmem')  # ~ 1G/core ~ 14.4G
    if depends_on:
        j.depends_on(*depends_on)

    j.declare_resource_group(output={'tar': '{root}.tar'})
    j.command(
        f"""set -e

    # We've seen some GenomicsDB performance regressions related to intervals, 
    # so we're going to pretend we only have a single interval
    # using the --merge-input-intervals arg. There's no data in between since 
    # we didn't run HaplotypeCaller over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    
    # The batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!

    (while true; do df -h; pwd; du -sh $(dirname {j.output['tar']}); free -m; sleep 300; done) &

    echo "Adding {len(sample_names_to_add)} samples: {', '.join(sample_names_to_add)}"
    {f'echo "Skipping adding {len(sample_names_to_skip)} samples that are already in the DB: '
     f'{", ".join(sample_names_to_skip)}"' if sample_names_to_skip else ''}

    gatk --java-options -Xms{java_mem}g \\
      GenomicsDBImport \\
      {genomicsdb_param} \\
      --batch-size 50 \\
      -L {interval} \\
      --sample-name-map {sample_map} \\
      --reader-threads {java_mem} \\
      --merge-input-intervals \\
      --consolidate

    df -h; pwd; du -sh $(dirname {j.output['tar']}); free -m
    """
    )
    return j, sample_names_will_be_in_db


def _add_genotype_gvcfs_job(
    b: hb.Batch,
    genomicsdb_path: str,
    reference: hb.ResourceGroup,
    dbsnp: str,
    overwrite: bool,
    number_of_samples: int,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    interval: Optional[hb.ResourceFile] = None,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Run joint-calling on all samples in a genomics database
    """
    job_name = 'Joint genotyping: GenotypeGVCFs'
    if interval_idx is not None:
        job_name += f' {interval_idx}/{number_of_intervals}'

    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.cpu(2)
    j.memory('standard')  # ~ 4G/core ~ 8G
    # 4G (fasta+fai+dict) + 4G per sample divided by the number of intervals
    j.storage(f'{4 + number_of_samples * 4 // number_of_intervals}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""
set -o pipefail
set -ex

cd $(dirname {j.output_vcf})

export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS
gsutil -q cp -r {genomicsdb_path} .

(while true; do df -h; pwd; free -m; sleep 300; done) &

df -h; pwd; free -m

gatk --java-options -Xms8g \\
  GenotypeGVCFs \\
  -R {reference.base} \\
  -O {j.output_vcf['vcf.gz']} \\
  -D {dbsnp} \\
  --only-output-calls-starting-in-intervals \\
  -V gendb://{basename(genomicsdb_path)} \\
  {f'-L {interval} ' if interval else ''} \\
  --merge-input-intervals

df -h; pwd; free -m
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_gnarly_genotyper_job(
    b: hb.Batch,
    genomicsdb_path: str,
    reference: hb.ResourceGroup,
    dbsnp: str,
    overwrite: bool,
    number_of_samples: int,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    interval: Optional[hb.ResourceGroup] = None,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Runs GATK GnarlyGenotyper on a combined_gvcf VCF bgzipped file.

    GnarlyGenotyper performs "quick and dirty" joint genotyping on large cohorts,
    pre-called with HaplotypeCaller, and post-processed with ReblockGVCF.

    HaplotypeCaller must be used with `-ERC GVCF` or `-ERC BP_RESOLUTION` to add
    genotype likelihoods.

    ReblockGVCF must be run to add all the annotations necessary for VQSR:
    QUALapprox, VarDP, RAW_MQandDP.

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: GnarlyGenotyper'
    if interval_idx is not None:
        job_name += f' {interval_idx}/{number_of_intervals}'

    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.cpu(2)
    j.memory('standard')  # ~ 4G/core ~ 8G
    # 4G (fasta+fai+dict) + 4G per sample divided by the number of intervals
    j.storage(f'{4 + number_of_samples * 4 // number_of_intervals}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    j.command(
        f"""
set -o pipefail
set -ex

cd $(dirname {j.output_vcf})

export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS
gsutil -q cp -r {genomicsdb_path} .
    
(while true; do df -h; pwd; free -m; sleep 300; done) &

df -h; pwd; free -m

gatk --java-options -Xms8g \\
  GnarlyGenotyper \\
  -R {reference.base} \\
  -O {j.output_vcf['vcf.gz']} \\
  -D {dbsnp} \\
  --only-output-calls-starting-in-intervals \\
  --keep-all-sites \\
  -V gendb://{basename(genomicsdb_path)} \\
  {f'-L {interval} ' if interval else ''} \\
  --create-output-variant-index

df -h; pwd; free -m
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_exccess_het_filter(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    excess_het_threshold: float = 54.69,
    interval: Optional[hb.ResourceGroup] = None,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Filter a large cohort callset on Excess Heterozygosity.

    The filter applies only to large callsets (`not is_small_callset`)

    Requires all samples to be unrelated.

    ExcessHet estimates the probability of the called samples exhibiting excess
    heterozygosity with respect to the null hypothesis that the samples are unrelated.
    The higher the score, the higher the chance that the variant is a technical artifact
    or that there is consanguinuity among the samples. In contrast to Inbreeding
    Coefficient, there is no minimal number of samples for this annotation.

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: ExcessHet filter'
    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'32G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail

    # Captring stderr to avoid Batch pod from crashing with OOM from millions of
    # warning messages from VariantFiltration, e.g.:
    # > JexlEngine - ![0,9]: 'ExcessHet > 54.69;' undefined variable ExcessHet
    gatk --java-options -Xms3g \\
      VariantFiltration \\
      --filter-expression 'ExcessHet > {excess_het_threshold}' \\
      --filter-name ExcessHet \\
      {f'-L {interval} ' if interval else ''} \\
      -O {j.output_vcf['vcf.gz']} \\
      -V {input_vcf['vcf.gz']} \\
      2> {j.stderr}
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_make_sites_only_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Create sites-only VCF with only site-level annotations.
    Speeds up the analysis in the AS-VQSR modeling step.

    Returns: a Job object with a single output j.sites_only_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: MakeSitesOnlyVcf'
    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'32G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      MakeSitesOnlyVcf \\
      -I {input_vcf['vcf.gz']} \\
      -O {j.output_vcf['vcf.gz']}
      """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_final_gather_vcf_job(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceGroup],
    overwrite: bool,
    output_vcf_path: str = None,
) -> Job:
    """
    Combines per-interval scattered VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    job_name = f'Gather {len(input_vcfs)} VCFs'
    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage(f'{1 + len(input_vcfs) * 1}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    (while true; do df -h; pwd free -m; sleep 300; done) &

    # --ignore-safety-checks makes a big performance difference so we include it in 
    # our invocation. This argument disables expensive checks that the file headers 
    # contain the same set of genotyped samples and that files are in order 
    # by position of first record.
    gatk --java-options -Xms{java_mem}g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      {input_cmdl} \\
      --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}
    
    df -h; pwd; free -m
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j
