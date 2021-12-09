#!/usr/bin/env python3

"""
Functions that craete Batch jobs that get from raw data to a combiner-ready GVCF
"""

import logging
from collections import defaultdict
from os.path import join, dirname, splitext, basename
from typing import Optional, List, Tuple

import pandas as pd
import hailtop.batch as hb
from hailtop.batch.job import Job
from joint_calling import sm_utils, resources
from joint_calling import utils

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def add_pre_combiner_jobs(
    b: hb.Batch,
    samples_df: pd.DataFrame,
    pre_combiner_bucket: str,
    output_suffix: str,
    overwrite: bool,
    analysis_project: str,  # pylint: disable=unused-argument
    force_samples: List[str],
    assume_gvcfs_are_ready: bool = False,
) -> Tuple[pd.DataFrame, str, List[Job]]:
    """
    Add jobs that prepare GVCFs for the combiner, if needed.

    Returns a Tuple of: a DataFrame with the sample metadata, a TSV file
    corresponding to that dataframe, and a list of jobs to wait for before
    submitting the combiner job
    """

    logger.info(f'Samples DF:\n{samples_df}')
    jobs_by_sample = defaultdict(list)

    def get_project_bucket(_proj):
        if _proj in ['syndip', 'giab']:
            proj_bucket = f'gs://cpg-reference/validation/{_proj}'
        else:
            proj_bucket = f'gs://cpg-{_proj}-{output_suffix}'
        return proj_bucket

    # Samples for which a raw GVCF is provided as input:
    gvcf_df = samples_df[samples_df.topostproc_gvcf != '-']
    for sn, proj, source, external_id, gvcf_path in zip(
        gvcf_df.s,
        gvcf_df.project,
        gvcf_df.source,
        gvcf_df.external_id,
        gvcf_df.topostproc_gvcf,
    ):
        proj_bucket = get_project_bucket(proj)
        gvcf_bucket = join(proj_bucket, 'gvcf')
        if source != '-':
            gvcf_bucket = join(gvcf_bucket, source)
        output_gvcf_path = join(gvcf_bucket, f'{sn}.g.vcf.gz')

        if sn in force_samples or (
            not assume_gvcfs_are_ready and 
            not utils.can_reuse(output_gvcf_path, overwrite)
        ):
            j = _add_postproc_gvcf_jobs(
                b=b,
                gvcf_path=gvcf_path,
                output_path=output_gvcf_path,
                sample_name=sn,
                project_name=proj,
                external_id=external_id,
            )
            jobs_by_sample[sn].append(j)
        samples_df.loc[sn, ['gvcf']] = output_gvcf_path
        logger.info(f'Updating sample {sn} gvcf to {output_gvcf_path}')

    # Samples for which a CRAM/BAM is provided as input for realignment:
    realign_df = samples_df[samples_df.realign_cram != '-']
    for sn, proj, source, input_cram, input_crai in zip(
        realign_df.s,
        realign_df.project,
        realign_df.source,
        realign_df.realign_cram,
        realign_df.realign_crai,
    ):
        assert isinstance(input_crai, str), (sn, input_crai)

        proj_bucket = get_project_bucket(proj)
        cram_bucket = join(proj_bucket, 'cram')
        if source != '-':
            cram_bucket = join(cram_bucket, source)
        output_cram_fpath = join(cram_bucket, f'{sn}.cram')
        output_crai_fpath = join(cram_bucket, f'{sn}.cram.crai')

        if not utils.can_reuse(output_cram_fpath, overwrite):
            alignment_input = sm_utils.AlignmentInput(
                bam_or_cram_path=input_cram,
                index_path=input_crai,
            )
            cram_j = _add_realign_jobs(
                b=b,
                output_path=output_cram_fpath,
                sample_name=sn,
                project_name=proj,
                alignment_input=alignment_input,
            )
            jobs_by_sample[sn].append(cram_j)
        samples_df.loc[sn, ['cram']] = output_cram_fpath
        samples_df.loc[sn, ['crai']] = output_crai_fpath

    # Samples for which a CRAM is provided as input for variant calling:
    hc_intervals_j = None
    cram_df = samples_df[samples_df.cram != '-']
    for sn, external_id, proj, source, input_cram, input_crai in zip(
        cram_df.s,
        cram_df.external_id,
        cram_df.project,
        cram_df.source,
        cram_df.cram,
        cram_df.crai,
    ):
        assert isinstance(input_crai, str), (sn, input_crai)

        proj_bucket = get_project_bucket(proj)
        gvcf_bucket = join(proj_bucket, 'gvcf')
        if source != '-':
            gvcf_bucket = join(gvcf_bucket, source)
        output_gvcf_path = join(gvcf_bucket, f'{sn}.g.vcf.gz')

        if not utils.can_reuse(output_gvcf_path, overwrite):
            if hc_intervals_j is None:
                hc_intervals_j = _add_split_intervals_job(
                    b=b,
                    interval_list=resources.UNPADDED_INTERVALS,
                    scatter_count=utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
                    ref_fasta=resources.REF_FASTA,
                )
            gvcf_j = _add_produce_gvcf_jobs(
                b=b,
                output_path=output_gvcf_path,
                sample_name=sn,
                project_name=proj,
                external_id=external_id,
                cram_fpath=input_cram,
                crai_fpath=input_crai,
                intervals_j=hc_intervals_j,
                tmp_bucket=join(pre_combiner_bucket, 'tmp'),
                overwrite=overwrite,
                depends_on=jobs_by_sample.get(sn, []),
            )
            jobs_by_sample[sn].append(gvcf_j)
        samples_df.loc[sn, ['gvcf']] = output_gvcf_path

    logger.info(f'Updated sample DF after post-processing GVCF inputs: {samples_df}')
    jobs = []
    for js in jobs_by_sample.values():
        jobs.extend(js)

    # Saving the resulting DataFrame as a TSV file
    gvcfs_tsv_path = join(pre_combiner_bucket, 'gvcfs.tsv')
    samples_df.to_csv(gvcfs_tsv_path, index=False, sep='\t', na_rep='NA')
    logger.info(f'Saved combiner-ready GVCF data to {gvcfs_tsv_path}')
    return samples_df, gvcfs_tsv_path, jobs


def _add_realign_jobs(
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
    j.cpu(total_cpu)
    j.memory('standard')
    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        }
    )

    reference = b.read_input_group(
        base=resources.REF_FASTA,
        fai=resources.REF_FASTA + '.fai',
        dict=resources.REF_FASTA.replace('.fasta', '')
        .replace('.fna', '')
        .replace('.fa', '')
        + '.dict',
        amb=resources.REF_FASTA + '.amb',
        ann=resources.REF_FASTA + '.ann',
        pac=resources.REF_FASTA + '.pac',
        o123=resources.REF_FASTA + '.0123',
        bwa2bit64=resources.REF_FASTA + '.bwt.2bit.64',
    )

    pull_inputs_cmd = ''

    if alignment_input.bam_or_cram_path:
        use_bazam = True
        bazam_cpu = 10
        bwa_cpu = 32
        bamsormadup_cpu = 10
        assert alignment_input.index_path
        assert not alignment_input.fqs1 and not alignment_input.fqs2
        j.storage(
            '400G' if alignment_input.bam_or_cram_path.endswith('.cram') else '1000G'
        )

        if alignment_input.bam_or_cram_path.startswith('gs://'):
            cram = b.read_input_group(
                base=alignment_input.bam_or_cram_path, index=alignment_input.index_path
            )
            cram_localized_path = cram.base
        else:
            # Can't use on Batch localization mechanism with `b.read_input_group`,
            # but have to manually localize with `wget`
            cram_name = basename(alignment_input.bam_or_cram_path)
            work_dir = dirname(j.output_cram.cram)
            cram_localized_path = join(work_dir, cram_name)
            index_ext = '.crai' if cram_name.endswith('.cram') else '.bai'
            crai_localized_path = join(work_dir, cram_name + index_ext)
            pull_inputs_cmd = (
                f'wget {alignment_input.bam_or_cram_path} -O {cram_localized_path}\n'
                f'wget {alignment_input.index_path} -O {crai_localized_path}'
            )
        r1_param = (
            f'<(bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base}'
            f' -n{bazam_cpu} -bam {cram_localized_path})'
        )
        r2_param = '-'
    else:
        assert alignment_input.fqs1 and alignment_input.fqs2
        use_bazam = False
        bwa_cpu = 32
        bamsormadup_cpu = 10
        j.storage('600G')

        files1 = [b.read_input(f1) for f1 in alignment_input.fqs1]
        files2 = [b.read_input(f1) for f1 in alignment_input.fqs2]
        r1_param = f'<(cat {" ".join(files1)})'
        r2_param = f'<(cat {" ".join(files2)})'
        logger.info(f'r1_param: {r1_param}')
        logger.info(f'r2_param: {r2_param}')

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

{pull_inputs_cmd}

bwa-mem2 mem -K 100000000 {'-p' if use_bazam else ''} -t{bwa_cpu} -Y \\
    -R '{rg_line}' {reference.base} {r1_param} {r2_param} | \\
bamsormadup inputformat=sam threads={bamsormadup_cpu} SO=coordinate \\
    M={j.duplicate_metrics} outputformat=sam \\
    tmpfile=$(dirname {j.output_cram.cram})/bamsormadup-tmp | \\
samtools view -T {reference.base} -O cram -o {j.output_cram.cram}

samtools index -@{total_cpu} {j.output_cram.cram} {j.output_cram['cram.crai']}

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


def _add_produce_gvcf_jobs(
    b: hb.Batch,
    output_path: str,
    sample_name: str,
    project_name: str,
    external_id: str,
    cram_fpath: str,
    crai_fpath: str,
    intervals_j: Job,
    tmp_bucket: str,
    overwrite: bool,
    depends_on: Optional[List[Job]] = None,
) -> Job:
    """
    Takes all samples with a 'file' of 'type'='bam' in `samples_df`,
    and runs HaplotypeCaller on them in parallel, then post-processes the GVCF

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
        base=resources.REF_FASTA,
        fai=resources.REF_FASTA + '.fai',
        dict=resources.REF_FASTA.replace('.fasta', '')
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
                cram_fpath=cram_fpath,
                crai_fpath=crai_fpath,
                interval=intervals_j.intervals[f'interval_{idx}'],
                reference=reference,
                interval_idx=idx,
                number_of_intervals=utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
                depends_on=depends_on,
            )
        )
    if depends_on:
        haplotype_caller_jobs[0].depends_on(*depends_on)

    hc_gvcf_path = join(tmp_bucket, 'haplotypecaller', f'{sample_name}.g.vcf.gz')
    merge_j = _add_merge_gvcfs_job(
        b=b,
        sample_name=sample_name,
        project_name=project_name,
        gvcfs=[j.output_gvcf for j in haplotype_caller_jobs],
        output_gvcf_path=hc_gvcf_path,
    )

    j = _add_postproc_gvcf_jobs(
        b=b,
        gvcf_path=hc_gvcf_path,
        output_path=output_path,
        sample_name=sample_name,
        project_name=project_name,
        external_id=external_id,
        depends_on=[merge_j],
    )
    return j


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
    cram_fpath: str,
    crai_fpath: str,
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
        job_name += f', {interval_idx}/{number_of_intervals}'
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

    export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json

    gatk --java-options "-Xms{java_mem}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
      HaplotypeCaller \\
      -R {reference.base} \\
      -I {cram_fpath} \\
      --read-index {crai_fpath} \\
      -L {interval} \\
      -O {j.output_gvcf['g.vcf.gz']} \\
      -G AS_StandardAnnotation \\
      -GQB 20 \
      -ERC GVCF \\

    df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m
    """
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
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


def _add_postproc_gvcf_jobs(
    b: hb.Batch,
    gvcf_path: str,
    output_path: str,
    sample_name: str,
    project_name: str,
    external_id: str,
    depends_on: Optional[List[Job]] = None,
) -> Job:
    """
    Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration.

    1. Subset GVCF to main chromosomes to avoid downstream errors
    2. Removes the DS INFO field that is added to some HGDP GVCFs to avoid errors
       from Hail about mismatched INFO annotations
    3. Renames sample name from external_sample_id to internal_sample_id
    """
    logger.info(f'Adding GVCF postproc job for sample {sample_name}, gvcf {gvcf_path}')

    j = b.new_job(f'{project_name}/{sample_name}: ReblockGVCF')
    j.image(utils.GATK_IMAGE)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'50G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )

    ref_fasta = resources.REF_FASTA
    ref_fai = resources.REF_FASTA + '.fai'
    ref_dict = (
        ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '') + '.dict'
    )

    j.command(f"""
    export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
    gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS

    function fail {{
      echo $1 >&2
      exit 1
    }}

    function retry {{
      local n=1
      local max=10
      local delay=30
      while true; do
        "$@" && break || {{
          if [[ $n -lt $max ]]; then
            ((n++))
            echo "Command failed. Attempt $n/$max:"
            sleep $delay;
          else
            fail "The command has failed after $n attempts."
          fi
        }}
      done
    }}

    # Retrying copying to avoid google bandwidth limits
    retry gsutil cp {gvcf_path} /io/batch/{sample_name}.g.vcf.gz
    retry gsutil cp {resources.NOALT_REGIONS} noalt-regions.bed

    # Copying reference data as well to avoid crazy logging costs 
    # for region requests
    retry gsutil cp {ref_fasta} /io/batch/{basename(ref_fasta)}
    retry gsutil cp {ref_fai}   /io/batch/{basename(ref_fai)}
    retry gsutil cp {ref_dict}  /io/batch/{basename(ref_dict)}

    # Reindexing just to make sure the index is not corrupted
    bcftools index --tbi /io/batch/{sample_name}.g.vcf.gz

    gatk --java-options "-Xms{mem_gb - 1}g" \\
    ReblockGVCF \\
    --reference /io/batch/{basename(ref_fasta)} \\
    -V /io/batch/{sample_name}.g.vcf.gz \\
    -do-qual-approx \\
    -O /io/batch/{sample_name}-reblocked.g.vcf.gz \\
    --create-output-variant-index true

    bcftools view /io/batch/{sample_name}-reblocked.g.vcf.gz -T noalt-regions.bed \\
    | bcftools annotate -x INFO/DS \\
    | bcftools reheader -s <(echo "{external_id} {sample_name}") \\
    | bcftools view -Oz -o {j.output_gvcf['g.vcf.gz']}

    bcftools index --tbi {j.output_gvcf['g.vcf.gz']}
    """)
    b.write_output(j.output_gvcf, output_path.replace('.g.vcf.gz', ''))
    if depends_on:
        j.depends_on(*depends_on)
    return j
