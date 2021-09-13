#!/usr/bin/env python3

"""
Functions that craete Batch jobs that get from raw data to a combiner-ready GVCF
"""

import logging
from os.path import join, dirname, splitext
from typing import Optional, List, Tuple

import pandas as pd
import hailtop.batch as hb
from hailtop.batch.job import Job
from joint_calling import sm_utils
from joint_calling import utils


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def add_pre_combiner_jobs(
    b: hb.Batch,
    samples_df: pd.DataFrame,
    pre_combiner_bucket: str,
    overwrite: bool,
    analysis_project: str,  # pylint: disable=unused-argument
) -> Tuple[pd.DataFrame, str, List[Job]]:
    """
    Add jobs that prepare GVCFs for the combiner, if needed.

    Returns a Tuple of: a DataFrame with the sample metadata, a TSV file
    corresponding to that dataframe, and a list of jobs to wait for before
    submitting the combiner job
    """

    gvcfs_tsv_path = join(pre_combiner_bucket, 'gvcfs.tsv')

    hc_intervals_j = None
    gvcf_jobs = []
    cram_df = samples_df[samples_df.cram.notna()]
    for s_id, proj, input_cram, input_crai in zip(
        cram_df.s, cram_df.project, cram_df.cram, cram_df.crai
    ):
        assert isinstance(input_crai, str)
        output_cram = join(pre_combiner_bucket, 'cram', f'{s_id}.cram')
        if not utils.can_reuse(output_cram, overwrite):
            alignment_input = sm_utils.AlignmentInput(
                bam_or_cram_path=input_cram,
                index_path=input_crai,
            )
            cram_j = _make_realign_jobs(
                b=b,
                output_path=output_cram,
                sample_name=s_id,
                project_name=proj,
                alignment_input=alignment_input,
            )
        else:
            cram_j = None
        samples_df.loc[s_id, ['cram']] = output_cram
        samples_df.loc[s_id, ['crai']] = output_cram + '.crai'

        output_gvcf_path = join(pre_combiner_bucket, 'gvcf', f'{s_id}.g.vcf.gz')
        if not utils.can_reuse(output_gvcf_path, overwrite):
            if hc_intervals_j is None:
                hc_intervals_j = _add_split_intervals_job(
                    b=b,
                    interval_list=utils.UNPADDED_INTERVALS,
                    scatter_count=utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
                    ref_fasta=utils.REF_FASTA,
                )
            gvcf_j = _make_produce_gvcf_jobs(
                b=b,
                output_path=output_gvcf_path,
                sample_name=s_id,
                project_name=proj,
                cram_path=output_cram,
                intervals_j=hc_intervals_j,
                tmp_bucket=join(pre_combiner_bucket, 'tmp'),
                overwrite=overwrite,
                depends_on=[cram_j] if cram_j else [],
            )
            gvcf_jobs.append(gvcf_j)
        samples_df.loc[s_id, ['gvcf']] = output_gvcf_path

    subset_gvcf_jobs, samples_df = _add_prep_gvcfs_for_combiner_steps(
        b=b,
        samples_df=samples_df,
        output_gvcf_bucket=join(pre_combiner_bucket, 'gvcf'),
        overwrite=overwrite,
        depends_on=gvcf_jobs,
    )
    samples_df.to_csv(gvcfs_tsv_path, index=False, sep='\t', na_rep='NA')
    logger.info(f'Saved combiner-ready GVCF data to {gvcfs_tsv_path}')
    return samples_df, gvcfs_tsv_path, subset_gvcf_jobs


def _make_realign_jobs(
    b: hb.Batch,
    output_path: str,
    sample_name: str,
    project_name: str,
    alignment_input: sm_utils.AlignmentInput,
) -> Job:
    """
    Runs BWA to realign reads back again to hg38.

    When the input is CRAM/BAM, uses Bazam to stream reads to BWA.
    """
    job_name = f'{project_name}/{sample_name}: BWA align'
    if utils.file_exists(output_path):
        return b.new_job(f'{job_name} [reuse]')

    logger.info(
        f'Parsing the "reads" metadata field and submitting the alignment '
        f'to write {output_path} for {sample_name}. '
    )
    j = b.new_job(job_name)
    j.image(utils.ALIGNMENT_IMAGE)
    total_cpu = 32

    reference = b.read_input_group(
        base=utils.REF_FASTA,
        fai=utils.REF_FASTA + '.fai',
        dict=utils.REF_FASTA.replace('.fasta', '')
        .replace('.fna', '')
        .replace('.fa', '')
        + '.dict',
        amb=utils.REF_FASTA + '.amb',
        ann=utils.REF_FASTA + '.ann',
        pac=utils.REF_FASTA + '.pac',
        o123=utils.REF_FASTA + '.0123',
        bwa2bit64=utils.REF_FASTA + '.bwt.2bit.64',
    )

    if alignment_input.bam_or_cram_path:
        use_bazam = True
        bazam_cpu = 10
        bwa_cpu = 32
        bamsormadup_cpu = 10

        assert alignment_input.index_path
        assert not alignment_input.fqs1 and not alignment_input.fqs2
        cram = b.read_input_group(
            base=alignment_input.bam_or_cram_path, index=alignment_input.index_path
        )
        r1_param = (
            f'<(bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base}'
            f' -n{bazam_cpu} -bam {cram.base})'
        )
        r2_param = '-'
    else:
        assert alignment_input.fqs1 and alignment_input.fqs2
        use_bazam = False
        bwa_cpu = 32
        bamsormadup_cpu = 10
        files1 = [b.read_input(f1) for f1 in alignment_input.fqs1]
        files2 = [b.read_input(f1) for f1 in alignment_input.fqs2]
        r1_param = f'<(cat {" ".join(files1)})'
        r2_param = f'<(cat {" ".join(files2)})'
        logger.info(f'r1_param: {r1_param}')
        logger.info(f'r2_param: {r2_param}')

    j.cpu(total_cpu)
    j.memory('standard')
    j.storage('750G')
    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'crai': '{root}.cram.crai',
        }
    )

    rg_line = f'@RG\\tID:{sample_name}\\tSM:{sample_name}'
    # BWA command options:
    # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
    # -p     smart pairing (ignoring in2.fq)
    # -t16   threads
    # -Y     use soft clipping for supplementary alignments
    # -R     read group header line such as '@RG\tID:foo\tSM:bar'
    command = f"""
set -o pipefail
set -ex

(while true; do df -h; pwd; du -sh $(dirname {j.output_cram.cram}); sleep 600; done) &

bwa-mem2 mem -K 100000000 {'-p' if use_bazam else ''} -t{bwa_cpu} -Y \\
    -R '{rg_line}' {reference.base} {r1_param} {r2_param} | \\
bamsormadup inputformat=sam threads={bamsormadup_cpu} SO=coordinate \\
    M={j.duplicate_metrics} outputformat=sam \\
    tmpfile=$(dirname {j.output_cram.cram})/bamsormadup-tmp | \\
samtools view -T {reference.base} -O cram -o {j.output_cram.cram}

samtools index -@{total_cpu} {j.output_cram.cram} {j.output_cram.crai}

df -h; pwd; du -sh $(dirname {j.output_cram.cram})
    """
    j.command(command)
    b.write_output(j.output_cram, splitext(output_path)[0])
    b.write_output(
        j.duplicate_metrics,
        join(
            dirname(output_path),
            'duplicate-metrics',
            f'{sample_name}-duplicate-metrics.csv',
        ),
    )
    return j


def _make_produce_gvcf_jobs(
    b: hb.Batch,
    output_path: str,
    sample_name: str,
    project_name: str,
    cram_path: str,
    intervals_j: Job,
    tmp_bucket: str,  # pylint: disable=unused-argument
    overwrite: bool,  # pylint: disable=unused-argument
    depends_on: Optional[List[Job]] = None,
) -> Job:
    """
    Takes all samples with a 'file' of 'type'='bam' in `samples_df`,
    and runs HaplotypeCaller on them, and sets a new 'file' of 'type'='gvcf'

    HaplotypeCaller is run in an interval-based sharded way, with per-interval
    HaplotypeCaller jobs defined in a nested loop.
    """
    job_name = f'{project_name}/{sample_name}: make GVCF'
    if utils.file_exists(output_path):
        return b.new_job(f'{job_name} [reuse]')
    logger.info(
        f'Submitting the variant calling jobs to write {output_path} for {sample_name}'
    )

    reference = b.read_input_group(
        base=utils.REF_FASTA,
        fai=utils.REF_FASTA + '.fai',
        dict=utils.REF_FASTA.replace('.fasta', '')
        .replace('.fna', '')
        .replace('.fa', '')
        + '.dict',
    )

    haplotype_caller_jobs = []
    for idx in range(utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS):
        haplotype_caller_jobs.append(
            _add_haplotype_caller_job(
                b,
                sample_name=sample_name,
                project_name=project_name,
                cram=b.read_input_group(
                    **{
                        'cram': cram_path,
                        'crai': cram_path + '.crai',
                    }
                ),
                interval=intervals_j.intervals[f'interval_{idx}'],
                reference=reference,
                interval_idx=idx,
                number_of_intervals=utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
                depends_on=depends_on,
            )
        )
    hc_gvcf_path = join(tmp_bucket, 'haplotypecaller', f'{sample_name}.g.vcf.gz')
    merge_j = _add_merge_gvcfs_job(
        b=b,
        sample_name=sample_name,
        project_name=project_name,
        gvcfs=[j.output_gvcf for j in haplotype_caller_jobs],
        output_gvcf_path=hc_gvcf_path,
    )

    if depends_on:
        haplotype_caller_jobs[0].depends_on(*depends_on)
    return merge_j


def _add_split_intervals_job(
    b: hb.Batch,
    interval_list: str,
    scatter_count: int,
    ref_fasta: str,
) -> Job:
    """
    Split genome into intervals to parallelise GnarlyGenotyper.

    Returns: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(f'Make {scatter_count} intervals')
    j.image(utils.GATK_IMAGE)
    java_mem = 3
    j.memory('standard')  # ~ 4G/core ~ 4G
    j.storage('16G')
    j.declare_resource_group(
        intervals={
            f'interval_{idx}': f'{{root}}/{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(scatter_count)
        }
    )

    j.command(
        f"""set -e

    # Modes other than INTERVAL_SUBDIVISION will produce an unpredicted number 
    # of intervals. But we have to expect exactly the {scatter_count} number of 
    # output files because our workflow is not dynamic.
    gatk --java-options -Xms{java_mem}g SplitIntervals \\
      -L {interval_list} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {ref_fasta} \\
      -mode INTERVAL_SUBDIVISION
      """
    )
    # Could save intervals to a bucket here to avoid rerunning the job
    return j


def _add_haplotype_caller_job(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    cram: hb.ResourceGroup,
    interval: hb.ResourceFile,
    reference: hb.ResourceGroup,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    depends_on: Optional[List[Job]] = None,
    output_gvcf_path: Optional[str] = None,
    overwrite: bool = False,
) -> Job:
    """
    Run HaplotypeCaller on an input BAM or CRAM, and output GVCF
    """
    job_name = f'{project_name}/{sample_name}: HaplotypeCaller'
    if interval_idx is not None:
        job_name += f', {sample_name} {interval_idx}/{number_of_intervals}'
    if utils.can_reuse(output_gvcf_path, overwrite):
        return b.new_job(f'{job_name} [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage('60G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )
    if depends_on:
        j.depends_on(*depends_on)

    j.command(
        f"""set -e
    (while true; do df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m; sleep 300; done) &

    gatk --java-options "-Xms{java_mem}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
      HaplotypeCaller \\
      -R {reference.base} \\
      -I {cram['cram']} \\
      -L {interval} \\
      -O {j.output_gvcf['g.vcf.gz']} \\
      -G AS_StandardAnnotation \\
      -GQB 20 \
      -ERC GVCF \\

    df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m
    """
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path)
    return j


def _add_merge_gvcfs_job(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    gvcfs: List[hb.ResourceGroup],
    output_gvcf_path: Optional[str],
) -> Job:
    """
    Combine by-interval GVCFs into a single sample GVCF file
    """

    job_name = f'{project_name}/{sample_name}: merge {len(gvcfs)} GVCFs'
    j = b.new_job(job_name)
    j.image(utils.PICARD_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage(f'{len(gvcfs) * 1.5 + 2}G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )

    input_cmd = ' '.join(f'INPUT={g["g.vcf.gz"]}' for g in gvcfs)

    j.command(
        f"""set -e

    (while true; do df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m; sleep 300; done) &

    java -Xms{java_mem}g -jar /usr/picard/picard.jar \
      MergeVcfs {input_cmd} OUTPUT={j.output_gvcf['g.vcf.gz']}

    df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m
      """
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def _add_prep_gvcfs_for_combiner_steps(
    b,
    samples_df: pd.DataFrame,
    output_gvcf_bucket: str,
    overwrite: bool,
    depends_on: Optional[List[Job]] = None,
) -> Tuple[List[Job], pd.DataFrame]:
    """
    Add steps required to prepare GVCFs from combining
    """
    noalt_regions = b.read_input(utils.NOALT_REGIONS)

    jobs = []
    logger.info(f'Samples DF: {samples_df}')
    for s_id, proj, external_id, gvcf_path in zip(
        samples_df.s, samples_df.project, samples_df.external_id, samples_df.gvcf
    ):
        output_gvcf_path = join(output_gvcf_bucket, f'{s_id}.g.vcf.gz')
        if not utils.can_reuse(output_gvcf_path, overwrite):
            logger.info(
                f'Adding reblock and subset jobs for sample {s_id}, gvcf {gvcf_path}'
            )
            gvcf = b.read_input_group(
                **{'g.vcf.gz': gvcf_path, 'g.vcf.gz.tbi': gvcf_path + '.tbi'}
            )
            reblock_j = _add_reblock_gvcf_job(
                b=b,
                sample_name=s_id,
                project_name=proj,
                input_gvcf=gvcf,
                overwrite=overwrite,
            )
            if depends_on:
                reblock_j.depends_on(*depends_on)
            jobs.append(
                _add_subset_noalt_step(
                    b=b,
                    sample_name=s_id,
                    project_name=proj,
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


#
# def _make_postproc_gvcf_jobs(
#     b: hb.Batch,
#     sample_name: str,
#     project_name: str,
#     input_gvcf_path: str,
#     out_gvcf_path: str,
#     noalt_regions: hb.ResourceFile,
#     overwrite: bool,  # pylint: disable=unused-argument
#     depends_on: Optional[List[Job]] = None,
# ) -> Job:
#     reblock_gvcf_job = _add_reblock_gvcf_job(
#         b,
#         sample_name=sample_name,
#         project_name=project_name,
#         input_gvcf=b.read_input_group(
#             **{'g.vcf.gz': input_gvcf_path, 'g.vcf.gz.tbi': input_gvcf_path + '.tbi'}
#         ),
#         overwrite=overwrite,
#     )
#     if depends_on:
#         reblock_gvcf_job.depends_on(*depends_on)
#     subset_to_noalt_job = _add_subset_noalt_step(
#         b,
#         sample_name=sample_name,
#         project_name=project_name,
#         input_gvcf=reblock_gvcf_job.output_gvcf,
#         noalt_regions=noalt_regions,
#         overwrite=overwrite,
#         output_gvcf_path=out_gvcf_path,
#     )
#     return subset_to_noalt_job


def _add_reblock_gvcf_job(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    input_gvcf: hb.ResourceGroup,
    overwrite: bool,
    output_gvcf_path: Optional[str] = None,
) -> Job:
    """
    Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration
    """
    job_name = f'{project_name}/{sample_name}: ReblockGVCF'
    if utils.can_reuse(output_gvcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
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

    j.command(
        f"""
    gatk --java-options "-Xms{mem_gb - 1}g" \\
        ReblockGVCF \\
        -V {input_gvcf['g.vcf.gz']} \\
        -do-qual-approx \\
        -O {j.output_gvcf['g.vcf.gz']} \\
        --create-output-variant-index true"""
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def _add_subset_noalt_step(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
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
    job_name = f'{project_name}/{sample_name}: SubsetToNoalt'
    j = b.new_job(job_name)
    j.image(utils.BCFTOOLS_IMAGE)
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
