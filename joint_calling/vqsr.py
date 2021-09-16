"""
Create jobs to create and apply a VQSR model
"""

import os
from os.path import join
from typing import List, Optional, Dict
import logging
import hailtop.batch as hb
from hailtop.batch.job import Job
from analysis_runner import dataproc

from joint_calling import utils
from joint_calling.utils import can_reuse

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


SNP_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.8,
    99.6,
    99.5,
    99.4,
    99.3,
    99.0,
    98.0,
    97.0,
    90.0,
]
SNP_RECALIBRATION_ANNOTATION_VALUES = [
    'AS_QD',
    'AS_MQRankSum',
    'AS_ReadPosRankSum',
    'AS_FS',
    'AS_SOR',
    'AS_MQ',
]
INDEL_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.5,
    99.0,
    97.0,
    96.0,
    95.0,
    94.0,
    93.5,
    93.0,
    92.0,
    91.0,
    90.0,
]
INDEL_RECALIBRATION_ANNOTATION_VALUES = [
    'AS_FS',
    'AS_SOR',
    'AS_ReadPosRankSum',
    'AS_MQRankSum',
    'AS_QD',
]


def make_vqsr_jobs(
    b: hb.Batch,
    combined_mt_path: str,
    hard_filter_ht_path: str,
    meta_ht_path: str,
    gvcf_count: int,
    work_bucket: str,
    web_bucket: str,
    depends_on: Optional[List[Job]],
    vqsr_params_d: Dict,
    scatter_count: int,
    output_vcf_path: str,
    overwrite: bool,
) -> Job:
    """
    Add jobs that perform the allele-specific VQSR variant QC

    :param b: Batch object to add jobs to
    :param combined_mt_path: path to a Matrix Table combined with the Hail VCF combiner
    :param hard_filter_ht_path: path to HT with samples that failed QC
    :param meta_ht_path: path to HT with sample QC metadata
    :param gvcf_count: number of input samples. Can't read from combined_mt_path as it
           might not be yet genereated the point of Batch job submission
    :param work_bucket: bucket for intermediate files
    :param web_bucket: bucket for plots and evaluation results (exposed via http)
    :param depends_on: job that the created jobs should only run after
    :param vqsr_params_d: parameters for VQSR
    :param scatter_count: number of shards to patition data for scattering
    :param output_vcf_path: path to write final recalibrated VCF to
    :param overwrite: whether to not reuse existing intermediate and output files
    :return: a final Job, and a path to the VCF with VQSR annotations
    """

    # Reference files. All options have defaults.
    unpadded_intervals_path = os.path.join(
        utils.GATK_REF_BUCKET, 'hg38.even.handcurated.20k.intervals'
    )
    dbsnp_vcf = os.path.join(
        utils.GATK_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf'
    )
    dbsnp_vcf_index = os.path.join(
        utils.GATK_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf.idx'
    )
    hapmap_resource_vcf = os.path.join(utils.GATK_REF_BUCKET, 'hapmap_3.3.hg38.vcf.gz')
    hapmap_resource_vcf_index = os.path.join(
        utils.GATK_REF_BUCKET, 'hapmap_3.3.hg38.vcf.gz.tbi'
    )
    omni_resource_vcf = os.path.join(utils.GATK_REF_BUCKET, '1000G_omni2.5.hg38.vcf.gz')
    omni_resource_vcf_index = os.path.join(
        utils.GATK_REF_BUCKET, '1000G_omni2.5.hg38.vcf.gz.tbi'
    )
    one_thousand_genomes_resource_vcf = os.path.join(
        utils.GATK_REF_BUCKET, '1000G_phase1.snps.high_confidence.hg38.vcf.gz'
    )
    one_thousand_genomes_resource_vcf_index = os.path.join(
        utils.GATK_REF_BUCKET, '1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi'
    )
    mills_resource_vcf = os.path.join(
        utils.GATK_REF_BUCKET, 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
    )
    mills_resource_vcf_index = os.path.join(
        utils.GATK_REF_BUCKET, 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'
    )
    axiom_poly_resource_vcf = os.path.join(
        utils.GATK_REF_BUCKET,
        'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz',
    )
    axiom_poly_resource_vcf_index = os.path.join(
        utils.GATK_REF_BUCKET,
        'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi',
    )
    ref_fasta = utils.REF_FASTA
    ref_fasta_index = utils.REF_FASTA + '.fai'
    ref_dict = (
        ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '') + '.dict'
    )
    ref_fasta = b.read_input_group(
        base=ref_fasta,
        dict=ref_dict,
        fai=ref_fasta_index,
    )
    dbsnp_vcf = b.read_input_group(base=dbsnp_vcf, index=dbsnp_vcf_index)
    hapmap_resource_vcf = b.read_input_group(
        base=hapmap_resource_vcf, index=hapmap_resource_vcf_index
    )
    omni_resource_vcf = b.read_input_group(
        base=omni_resource_vcf, index=omni_resource_vcf_index
    )
    one_thousand_genomes_resource_vcf = b.read_input_group(
        base=one_thousand_genomes_resource_vcf,
        index=one_thousand_genomes_resource_vcf_index,
    )
    mills_resource_vcf = b.read_input_group(
        base=mills_resource_vcf, index=mills_resource_vcf_index
    )
    axiom_poly_resource_vcf = b.read_input_group(
        base=axiom_poly_resource_vcf, index=axiom_poly_resource_vcf_index
    )
    dbsnp_resource_vcf = dbsnp_vcf

    is_small_callset = gvcf_count < 1000
    # 1. For small callsets, we don't apply the ExcessHet filtering.
    # 2. For small callsets, we gather the VCF shards and collect QC metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics
    # collected from them.
    is_huge_callset = gvcf_count >= 100000
    # For huge callsets, we allocate more memory for the SNPs Create Model step

    small_disk = 30 if is_small_callset else (50 if not is_huge_callset else 100)
    medium_disk = 50 if is_small_callset else (100 if not is_huge_callset else 200)
    huge_disk = 100 if is_small_callset else (500 if not is_huge_callset else 2000)

    job_name = 'AS-VQSR: MT to VCF'
    combined_vcf_path = join(work_bucket, 'input.vcf.gz')
    if not can_reuse(combined_vcf_path, overwrite):
        mt_to_vcf_job = dataproc.hail_dataproc_job(
            b,
            f'{utils.SCRIPTS_DIR}/mt_to_vcf.py --overwrite '
            f'--mt {combined_mt_path} '
            f'--meta-ht {meta_ht_path} '
            f'--hard-filtered-samples-ht {hard_filter_ht_path} '
            f'-o {combined_vcf_path} ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=depends_on,
            # hl.export_vcf() uses non-preemptible workers' disk to merge VCF files.
            # 10 samples take 2.3G, 400 samples take 60G, which roughly matches
            # `huge_disk` (also used in the AS-VQSR VCF-gather job)
            worker_boot_disk_size=huge_disk,
            job_name=job_name,
        )
    else:
        mt_to_vcf_job = b.new_job(f'{job_name} [reuse]')

    split_intervals_job = add_split_intervals_step(
        b,
        unpadded_intervals_path,
        scatter_count,
        ref_fasta,
        disk_size=small_disk,
    )
    intervals = split_intervals_job.intervals

    tabix_job = add_tabix_step(b, combined_vcf_path, medium_disk)
    tabix_job.depends_on(mt_to_vcf_job)

    gathered_vcf = tabix_job.combined_vcf
    scattered_vcfs = [gathered_vcf for _ in range(scatter_count)]

    if not is_small_callset:
        # ExcessHet filtering applies only to callsets with a large number of samples,
        # e.g. hundreds of unrelated samples. Small cohorts should not trigger ExcessHet
        # filtering as values should remain small. Note cohorts of consanguinous samples
        # will inflate ExcessHet, and it is possible to limit the annotation to founders
        # for such cohorts by providing a pedigree file during variant calling.
        hard_filtered_vcfs = [
            add_hard_filter_step(
                b,
                input_vcf=gathered_vcf,
                interval=intervals[f'interval_{idx}'],
                excess_het_threshold=vqsr_params_d['min_excess_het'],
                disk_size=medium_disk,
            ).output_vcf
            for idx in range(scatter_count)
        ]
        scattered_vcfs = hard_filtered_vcfs

        gathered_vcf = add_sites_only_gather_vcf_step(
            b,
            input_vcfs=scattered_vcfs,
            disk_size=medium_disk,
        ).output_vcf

    indels_variant_recalibrator_job = add_indels_variant_recalibrator_step(
        b,
        sites_only_variant_filtered_vcf=gathered_vcf,
        mills_resource_vcf=mills_resource_vcf,
        axiom_poly_resource_vcf=axiom_poly_resource_vcf,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        disk_size=small_disk,
        web_bucket=web_bucket,
        work_bucket=web_bucket,
    )
    indels_recalibration = indels_variant_recalibrator_job.recalibration
    indels_tranches = indels_variant_recalibrator_job.tranches

    snp_max_gaussians = 6
    if is_small_callset:
        snp_max_gaussians = 4
    elif is_huge_callset:
        snp_max_gaussians = 8

    if is_huge_callset:
        # Run SNP recalibrator in a scattered mode
        model_file = add_snps_variant_recalibrator_create_model_step(
            b,
            sites_only_variant_filtered_vcf=gathered_vcf,
            hapmap_resource_vcf=hapmap_resource_vcf,
            omni_resource_vcf=omni_resource_vcf,
            one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
            dbsnp_resource_vcf=dbsnp_resource_vcf,
            disk_size=small_disk,
            web_bucket=web_bucket,
            work_bucket=work_bucket,
            is_small_callset=is_small_callset,
            is_huge_callset=is_huge_callset,
            max_gaussians=snp_max_gaussians,
        ).model_file
        # model_file = b.read_input('gs://playground-au/batch/859e9a/18/model_report')

        snps_recalibrator_jobs = [
            add_snps_variant_recalibrator_scattered_step(
                b,
                sites_only_vcf=scattered_vcfs[idx],
                interval=intervals[f'interval_{idx}'],
                model_file=model_file,
                hapmap_resource_vcf=hapmap_resource_vcf,
                omni_resource_vcf=omni_resource_vcf,
                one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
                dbsnp_resource_vcf=dbsnp_resource_vcf,
                disk_size=small_disk,
                max_gaussians=snp_max_gaussians,
            )
            for idx in range(scatter_count)
        ]
        snps_recalibrations = [j.recalibration for j in snps_recalibrator_jobs]
        snps_tranches = [j.tranches for j in snps_recalibrator_jobs]
        snps_gathered_tranches = add_snps_gather_tranches_step(
            b,
            tranches=snps_tranches,
            disk_size=small_disk,
        ).out_tranches

        scattered_vcfs = [
            add_apply_recalibration_step(
                b,
                input_vcf=scattered_vcfs[idx],
                interval=intervals[f'interval_{idx}'],
                indels_recalibration=indels_recalibration,
                indels_tranches=indels_tranches,
                snps_recalibration=snps_recalibrations[idx],
                snps_tranches=snps_gathered_tranches,
                disk_size=huge_disk,
                indel_filter_level=vqsr_params_d['indel_filter_level'],
                snp_filter_level=vqsr_params_d['snp_filter_level'],
            ).recalibrated_vcf
            for idx in range(scatter_count)
        ]
        recalibrated_gathered_vcf_job = _add_final_gather_vcf_step(
            b,
            input_vcfs=scattered_vcfs,
            disk_size=huge_disk,
        )
        recalibrated_gathered_vcf = recalibrated_gathered_vcf_job.output_vcf

    else:
        snps_recalibrator_job = add_snps_variant_recalibrator_step(
            b,
            sites_only_variant_filtered_vcf=gathered_vcf,
            hapmap_resource_vcf=hapmap_resource_vcf,
            omni_resource_vcf=omni_resource_vcf,
            one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
            dbsnp_resource_vcf=dbsnp_resource_vcf,
            disk_size=small_disk,
            max_gaussians=snp_max_gaussians,
            web_bucket=web_bucket,
            work_bucket=work_bucket,
        )
        snps_recalibration = snps_recalibrator_job.recalibration
        snps_tranches = snps_recalibrator_job.tranches

        recalibrated_gathered_vcf_job = add_apply_recalibration_step(
            b,
            input_vcf=gathered_vcf,
            indels_recalibration=indels_recalibration,
            indels_tranches=indels_tranches,
            snps_recalibration=snps_recalibration,
            snps_tranches=snps_tranches,
            disk_size=huge_disk,
            indel_filter_level=vqsr_params_d['indel_filter_level'],
            snp_filter_level=vqsr_params_d['snp_filter_level'],
        )
        recalibrated_gathered_vcf = recalibrated_gathered_vcf_job.recalibrated_vcf

    final_gathered_vcf_job = _add_filter_sb_step(
        b,
        input_vcf=recalibrated_gathered_vcf,
        disk_size=medium_disk,
        output_vcf_path=output_vcf_path,
    )

    _add_variant_eval_step(
        b,
        input_vcf=final_gathered_vcf_job.output_vcf,
        ref_fasta=ref_fasta,
        dbsnp_vcf=dbsnp_vcf,
        output_path=os.path.join(web_bucket, 'variant-eval.txt'),
        disk_size=huge_disk,
    )

    return final_gathered_vcf_job


def add_tabix_step(
    b: hb.Batch,
    vcf_path: str,
    disk_size: int,
) -> Job:
    """
    Regzip and tabix the combined VCF (for some reason the one output with mt2vcf
    is not block-gzipped)
    """
    j = b.new_job('AS-VQSR: Tabix')
    j.image(utils.BCFTOOLS_IMAGE)
    j.memory(f'8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        combined_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    vcf_inp = b.read_input(vcf_path)
    j.command(
        f"""set -e
        gunzip {vcf_inp} -c | bgzip -c > {j.combined_vcf['vcf.gz']}
        tabix -p vcf {j.combined_vcf['vcf.gz']}
        """
    )
    return j


def add_split_intervals_step(
    b: hb.Batch,
    interval_list: hb.ResourceFile,
    scatter_count: int,
    ref_fasta: hb.ResourceGroup,
    disk_size: int,
) -> Job:
    """
    Split genome into intervals to parallelise GnarlyGenotyper.

    Returns: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job('AS-VQSR: SplitIntervals')
    j.image(utils.GATK_IMAGE)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        intervals={
            f'interval_{idx}': f'{{root}}/{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(scatter_count)
        }
    )

    j.command(
        f"""set -e

    # Modes other than INTERVAL_SUBDIVISION will produce an unpredictable number
    # of intervals. But we have to produce exactly {scatter_count} number of
    # output files because our workflow is not dynamic.
    gatk --java-options -Xms{mem_gb - 1}g SplitIntervals \\
      -L {interval_list} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {ref_fasta.base} \\
      -mode INTERVAL_SUBDIVISION
      """
    )
    return j


def add_gnarly_genotyper_on_vcf_step(
    b: hb.Batch,
    combined_gvcf: hb.ResourceGroup,
    ref_fasta: hb.ResourceGroup,
    dbsnp_vcf: hb.ResourceGroup,
    disk_size: int,
    interval: Optional[hb.ResourceGroup] = None,
) -> Job:
    """
    Runs GATK GnarlyGenotyper on a combined_gvcf VCF bgzipped file.

    GnarlyGenotyper performs "quick and dirty" joint genotyping on large cohorts,
    pre-called with HaplotypeCaller, and post-processed with ReblockGVCF.

    HaplotypeCaller must be used with `-ERC GVCF` or `-ERC BP_RESOLUTION` to add
    genotype likelihoods.

    ReblockGVCF must be run to remove low quality variants, as well as to add all the
    annotations necessary for VQSR: QUALapprox, VarDP, RAW_MQandDP.

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    j = b.new_job('AS-VQSR: GnarlyGenotyperOnVcf')
    # GnarlyGenotyper crashes with NullPointerException when using standard GATK docker
    j.image(utils.GNARLY_IMAGE)
    j.memory(f'32G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -e

    gatk --java-options -Xms8g \\
      GnarlyGenotyper \\
      -R {ref_fasta.base} \\
      -O {j.output_vcf['vcf.gz']} \\
      -D {dbsnp_vcf.base} \\
      --only-output-calls-starting-in-intervals \\
      --keep-all-sites \\
      -V {combined_gvcf['vcf.gz']} \\
      {f'-L {interval} ' if interval else ''} \\
      --create-output-variant-index"""
    )
    return j


def add_hard_filter_step(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    excess_het_threshold: float,
    disk_size: int,
    interval: Optional[hb.ResourceGroup] = None,
) -> Job:
    """
    Hard-filter a large cohort callset on Excess Heterozygosity.

    Applies only to large callsets (`not is_small_callset`)

    Requires all samples to be unrelated.

    ExcessHet estimates the probability of the called samples exhibiting excess
    heterozygosity with respect to the null hypothesis that the samples are unrelated.
    The higher the score, the higher the chance that the variant is a technical artifact
    or that there is consanguinuity among the samples. In contrast to Inbreeding
    Coefficient, there is no minimal number of samples for this annotation.

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    j = b.new_job('AS-VQSR: HardFilter')
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'{disk_size}G')
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
    return j


def add_make_sites_only_vcf_step(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    disk_size: int,
    interval: Optional[hb.ResourceGroup] = None,
) -> Job:
    """
    Create sites-only VCF with only site-level annotations.
    Speeds up the analysis in the modeling step.

    Returns: a Job object with a single output j.sites_only_vcf of type ResourceGroup
    """
    j = b.new_job('AS-VQSR: MakeSitesOnlyVcf')
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        sites_only_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      MakeSitesOnlyVcf \\
      -I {input_vcf['vcf.gz']} \\
      -O {j.sites_only_vcf['vcf.gz']} \\
      {f'-L {interval} ' if interval else ''}"""
    )
    return j


def add_sites_only_gather_vcf_step(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceFile],
    disk_size: int,
) -> Job:
    """
    Gathers VCF files from scattered operations into a single VCF file

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    j = b.new_job('AS-VQSR: SitesOnlyGatherVcf')
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'{disk_size}G')

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in
    # our invocation. This argument disables expensive checks that the file headers
    # contain the same set of genotyped samples and that files are in order by position
    # of first record.
    gatk --java-options -Xms6g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      {input_cmdl} \\
      --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}"""
    )
    return j


def add_indels_variant_recalibrator_step(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    mills_resource_vcf: hb.ResourceGroup,
    axiom_poly_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    web_bucket: str = None,
    work_bucket: str = None,
    max_gaussians: int = 4,
) -> Job:
    """
    Run VariantRecalibrator to calculate VQSLOD tranches for indels

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value. 4 is a
    reasonable default for indels, as their number is smaller than SNPs.

    Returns: a Job object with 3 outputs: j.recalibration (ResourceGroup), j.tranches,
    and j.indel_rscript_file. The latter is usedful to produce the optional tranche plot.
    """
    j = b.new_job('AS-VQSR: IndelsVariantRecalibrator')
    j.image(utils.GATK_IMAGE)
    mem_gb = 32
    j.memory(f'{mem_gb}G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join(
        [f'-tranche {v}' for v in INDEL_RECALIBRATION_TRANCHE_VALUES]
    )
    an_cmdl = ' '.join([f'-an {v}' for v in INDEL_RECALIBRATION_ANNOTATION_VALUES])
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{mem_gb - 1}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode INDEL \\
      --use-allele-specific-annotations \\
      --max-gaussians {max_gaussians} \\
      -resource:mills,known=false,training=true,truth=true,prior=12 {mills_resource_vcf.base} \\
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {axiom_poly_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 {dbsnp_resource_vcf.base} \\
      --rscript-file {j.indel_rscript_file}
      
      ls $(dirname {j.indel_rscript_file})

      ln {j.indel_rscript_file}.pdf {j.indel_features_pdf}
      """
    )
    if work_bucket:
        b.write_output(
            j.indel_rscript_file,
            os.path.join(work_bucket, 'recalibration-indels-features.Rscript'),
        )
    if web_bucket:
        b.write_output(
            j.indel_features_pdf,
            os.path.join(web_bucket, 'recalibration-indels-features.pdf'),
        )
    return j


def add_snps_variant_recalibrator_create_model_step(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    web_bucket: str = None,
    work_bucket: str = None,
    is_small_callset: bool = False,
    is_huge_callset: bool = False,
    max_gaussians: int = 4,
) -> Job:
    """
    First step of VQSR for SNPs: run VariantRecalibrator to subsample variants
    and produce a file of the VQSR model.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    Returns: a Job object with 2 outputs: j.model and j.snp_rscript_file.
    The latter is useful to produce the optional tranche plot.
    """
    j = b.new_job('AS-VQSR: SNPsVariantRecalibratorCreateModel')
    j.image(utils.GATK_IMAGE)
    mem_gb = 64 if not is_small_callset else 128
    j.memory(f'{mem_gb}G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    downsample_factor = 75 if is_huge_callset else 10

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join([f'-an {v}' for v in SNP_RECALIBRATION_ANNOTATION_VALUES])
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{mem_gb - 2}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      --use-allele-specific-annotations \\
      --sample-every-Nth-variant {downsample_factor} \\
      --output-model {j.model_file} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base} \\
      --rscript-file {j.snp_rscript}

      ls $(dirname {j.snp_rscript})

      ln {j.snp_rscript}.pdf {j.snp_rscript_pdf}
      ln {j.tranches}.pdf {j.tranches_pdf}
      """
    )
    if work_bucket:
        b.write_output(
            j.snp_rscript,
            os.path.join(work_bucket, 'recalibration-snps-features.RScript'),
        )
    if web_bucket:
        b.write_output(
            j.snp_rscript_pdf,
            os.path.join(web_bucket, 'recalibration-snps-features.pdf'),
        )
        b.write_output(
            j.tranches_pdf,
            os.path.join(web_bucket, 'recalibration-snps-tranches.pdf'),
        )
    return j


def add_snps_variant_recalibrator_scattered_step(
    b: hb.Batch,
    sites_only_vcf: hb.ResourceGroup,
    model_file: hb.ResourceGroup,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    interval: Optional[hb.ResourceGroup] = None,
    max_gaussians: int = 4,
) -> Job:
    """
    Second step of VQSR for SNPs: run VariantRecalibrator scattered to apply
    the VQSR model file to each genomic interval.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    Returns: a Job object with 2 outputs: j.recalibration (ResourceGroup) and j.tranches
    """
    j = b.new_job('AS-VQSR: SNPsVariantRecalibratorScattered')

    j.image(utils.GATK_IMAGE)
    mem_gb = 64  # ~ twice the sum of all input resources and input VCF sizes
    j.memory(f'{mem_gb}G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join([f'-an {v}' for v in SNP_RECALIBRATION_ANNOTATION_VALUES])
    j.command(
        f"""set -euo pipefail

    MODEL_REPORT={model_file}

    gatk --java-options -Xms{mem_gb - 1}g \\
      VariantRecalibrator \\
      -V {sites_only_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      {f'-L {interval} ' if interval else ''} \\
      --use-allele-specific-annotations \\
      --input-model {model_file} --output-tranches-for-scatter \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base}"""
    )
    return j


def add_snps_variant_recalibrator_step(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    web_bucket: str,
    work_bucket: str,
    disk_size: int,
    max_gaussians: int = 4,
) -> Job:
    """
    Recalibrate SNPs in one run (alternative to scatter-gather approach)
    """
    j = b.new_job('AS-VQSR: SNPsVariantRecalibrator')

    j.image(utils.GATK_IMAGE)
    mem_gb = 64  # ~ twice the sum of all input resources and input VCF sizes
    j.memory(f'{mem_gb}G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join([f'-an {v}' for v in SNP_RECALIBRATION_ANNOTATION_VALUES])
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{mem_gb - 1}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      --use-allele-specific-annotations \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base} \\
      --rscript-file {j.snp_rscript}

      ln {j.snp_rscript}.pdf {j.snp_rscript_pdf}
      ln {j.tranches}.pdf {j.tranches_pdf}
      """
    )

    if work_bucket:
        b.write_output(
            j.snp_rscript,
            os.path.join(work_bucket, 'recalibration-snps-features.RScript'),
        )
    if web_bucket:
        b.write_output(
            j.snp_rscript_pdf,
            os.path.join(web_bucket, 'recalibration-snps-features.pdf'),
        )
        b.write_output(
            j.tranches_pdf,
            os.path.join(web_bucket, 'recalibration-snps-tranches.pdf'),
        )
    return j


def add_snps_gather_tranches_step(
    b: hb.Batch,
    tranches: List[hb.ResourceFile],
    disk_size: int,
) -> Job:
    """
    Third step of VQSR for SNPs: run GatherTranches to gather scattered per-interval
    tranches outputs.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches

    Returns: a Job object with one output j.out_tranches
    """
    j = b.new_job('AS-VQSR: SNPGatherTranches')
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    inputs_cmdl = ' '.join([f'--input {t}' for t in tranches])
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      GatherTranches \\
      --mode SNP \\
      {inputs_cmdl} \\
      --output {j.out_tranches}"""
    )
    return j


def add_apply_recalibration_step(
    b: hb.Batch,
    input_vcf: hb.ResourceFile,
    indels_recalibration: hb.ResourceGroup,
    indels_tranches: hb.ResourceFile,
    snps_recalibration: hb.ResourceGroup,
    snps_tranches: hb.ResourceFile,
    disk_size: int,
    indel_filter_level: float,
    snp_filter_level: float,
    interval: Optional[hb.ResourceGroup] = None,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Apply a score cutoff to filter variants based on a recalibration table.
    Runs ApplyVQSR twice to apply first indel, then SNP recalibrations.

    Targets indel_filter_level and snp_filter_level sensitivities. The tool matches
    them internally to a VQSLOD score cutoff based on the model's estimated sensitivity
    to a set of true variants.

    The filter determination is not just a pass/fail process. The tool evaluates for
    each variant which "tranche", or slice of the dataset, it falls into in terms of
    sensitivity to the truthset. Variants in tranches that fall below the specified
    truth sensitivity filter level have their FILTER field annotated with the
    corresponding tranche level. This results in a callset that is filtered to the
    desired level but retains the information necessary to increase sensitivity
    if needed.

    Returns: a Job object with one ResourceGroup output j.recalibrated_vcf, correponding
    to a VCF with tranche annotated in the FILTER field
    """
    j = b.new_job('AS-VQSR: ApplyRecalibration')
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        recalibrated_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail
        
    df -h; pwd; du -sh $(dirname {j.recalibrated_vcf['vcf.gz']})
    
    TMP_DIR=$(dirname {j.recalibrated_vcf['vcf.gz']})/tmp
    mkdir $TMP_DIR

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      --tmp-dir $TMP_DIR \\
      -O tmp.indel.recalibrated.vcf.gz \\
      -V {input_vcf['vcf.gz']} \\
      --recal-file {indels_recalibration} \\
      --tranches-file {indels_tranches} \\
      --truth-sensitivity-filter-level {indel_filter_level} \\
      --create-output-variant-index true \\
      -mode INDEL \\
      {f'-L {interval} ' if interval else ''} \\
      --use-allele-specific-annotations

    df -h; pwd; du -sh $(dirname {j.recalibrated_vcf['vcf.gz']})

    rm {input_vcf['vcf.gz']} {indels_recalibration} {indels_tranches}
    rm -rf $TMP_DIR
    mkdir $TMP_DIR

    df -h; pwd; du -sh $(dirname {j.recalibrated_vcf['vcf.gz']})

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      --tmp-dir $TMP_DIR \\
      -O {j.recalibrated_vcf['vcf.gz']} \\
      -V tmp.indel.recalibrated.vcf.gz \\
      --recal-file {snps_recalibration} \\
      --tranches-file {snps_tranches} \\
      --truth-sensitivity-filter-level {snp_filter_level} \\
      --create-output-variant-index true \\
      -mode SNP \\
      {f'-L {interval} ' if interval else ''} \\
      --use-allele-specific-annotations

    df -h; pwd; du -sh $(dirname {j.recalibrated_vcf['vcf.gz']})
    """
    )

    if output_vcf_path:
        b.write_output(j.recalibrated_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j


def add_collect_metrics_sharded_step(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    dbsnp_vcf: hb.ResourceGroup,
    interval_list: hb.ResourceFile,
    ref_dict: hb.ResourceFile,
    disk_size: int,
):
    """
    Run CollectVariantCallingMetrics for site-level evaluation.

    This produces detailed and summary metrics report files. The summary metrics
    provide cohort-level variant metrics and the detailed metrics segment variant
    metrics for each sample in the callset. The detail metrics give the same metrics
    as the summary metrics for the samples plus several additional metrics.

    These are explained in detail at
    https://broadinstitute.github.io/picard/picard-metric-definitions.html.

    Returns: Job object with a single ResourceGroup output j.metrics, with
    j.metrics.detail_metrics and j.metrics.summary_metrics ResourceFiles
    """
    j = b.new_job('AS-VQSR: CollectMetricsSharded')
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.cpu(2)
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        metrics={
            'detail_metrics': '{root}.variant_calling_detail_metrics',
            'summary_metrics': '{root}.variant_calling_summary_metrics',
        }
    )

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      CollectVariantCallingMetrics \\
      --INPUT {input_vcf['vcf.gz']} \\
      --DBSNP {dbsnp_vcf.base} \\
      --SEQUENCE_DICTIONARY {ref_dict} \\
      --OUTPUT {j.metrics} \\
      --THREAD_COUNT 8 \\
      --TARGET_INTERVALS {interval_list}"""
    )
    return j


def _add_final_gather_vcf_step(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceGroup],
    disk_size: int,
    output_vcf_path: str = None,
) -> Job:
    """
    Combines recalibrated VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    j = b.new_job('AS-VQSR: FinalGatherVcf')
    j.image(utils.GATK_IMAGE)
    j.memory(f'8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in 
    # our invocation. This argument disables expensive checks that the file headers 
    # contain the same set of genotyped samples and that files are in order 
    # by position of first record.
    gatk --java-options -Xms6g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      {input_cmdl} \\
      --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}"""
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j


def _add_filter_sb_step(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    disk_size: int,
    output_vcf_path: str = None,
) -> Job:
    """
    Removes the INFO/SB field from a VCF.

    The reason we are doing that is because gatk ApplyVQSR replaces the VCF header
    ##INFO=<ID=SB,Number=.,Type=Int,Description="Strand Bias">

    with
    ##INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">

    Even though the actual SB field is still a list of integers: SB=5,2,18,29
    It breaks parsing the VCF into Hail.
    """
    j = b.new_job('AS-VQSR: Remove SB')
    j.image(utils.BCFTOOLS_IMAGE)
    j.memory(f'8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""
    bcftools annotate -x INFO/SB {input_vcf['vcf.gz']} -Oz -o {j.output_vcf['vcf.gz']} && tabix {j.output_vcf['vcf.gz']}
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j


def _add_variant_eval_step(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    ref_fasta: hb.ResourceGroup,
    dbsnp_vcf: hb.ResourceGroup,
    disk_size: int,
    output_path: str = None,
) -> Job:
    """
    Run VariantEval for site-level evaluation.
    Saves the QC to `output_path` bucket
    """
    j = b.new_job('AS-VQSR: VariantEval')
    j.image(utils.GATK_IMAGE)
    j.memory(f'8G')
    j.storage(f'{disk_size}G')

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      VariantEval \\
      --eval {input_vcf['vcf.gz']} \\
      -R {ref_fasta.base} \\
      -D {dbsnp_vcf.base} \\
      --output {j.output}"""
    )
    if output_path:
        b.write_output(j.output, output_path)
    return j


def add_gather_variant_calling_metrics_step(
    b: hb.Batch,
    input_details: List[hb.ResourceGroup],
    input_summaries: List[hb.ResourceGroup],
    disk_size: int,
    output_path_prefix: str = None,
) -> Job:
    """
    Combines metrics from multiple CollectVariantCallingMetrics runs.

    Returns: Job object with a single ResourceGroup output j.metrics, with
    j.metrics.detail_metrics and j.metrics.summary_metrics ResourceFiles

    Saves the QC results to a bucket with the `output_path_prefix` prefix
    """
    j = b.new_job('AS-VQSR: GatherVariantCallingMetrics')
    j.image(utils.GATK_IMAGE)
    j.memory(f'8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        metrics={
            'detail_metrics': '{root}.variant_calling_detail_metrics',
            'summary_metrics': '{root}.variant_calling_summary_metrics',
        }
    )

    input_cmdl = ' '.join('--INPUT {f} ' for f in input_details + input_summaries)
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms2g \\
      AccumulateVariantCallingMetrics \\
      {input_cmdl} \\
      --OUTPUT {j.metrics}"""
    )
    if output_path_prefix:
        b.write_output(j.metrics, output_path_prefix)
    return j
