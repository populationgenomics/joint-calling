import math, os
from typing import Union, Optional, List

import click
import hailtop.batch as hb


GATK_VERSION = "4.2.0.0"
GATK_DOCKER = f"us.gcr.io/broad-gatk/gatk:{GATK_VERSION}"
# GnarlyGenotyper crashes with NullPointerException when using GATK docker
GNARLY_DOCKER = "gcr.io/broad-dsde-methods/gnarly_genotyper:hail_ukbb_300K"


@click.command()
@click.option("--combined_gvcf", "combined_gvcf", type=str, required=True)
@click.option("--output_bucket", "output_bucket", type=str, required=True)
@click.option("--num_gvcfs", "num_gvcfs", type=int, required=True)
@click.option("--callset_name", "callset_name", type=str, required=True)
@click.option(
    "--unpadded_intervals_file",
    "unpadded_intervals_file",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/hg38.even.handcurated.20k.intervals",
)
@click.option(
    "--ref_fasta",
    "ref_fasta",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
)
@click.option(
    "--ref_fasta_index",
    "ref_fasta_index",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
)
@click.option(
    "--ref_dict",
    "ref_dict",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
)
@click.option(
    "--dbsnp_vcf",
    "dbsnp_vcf",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
)
@click.option(
    "--dbsnp_vcf_index",
    "dbsnp_vcf_index",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
)
@click.option("--small_disk", "small_disk", type=int, default=50)
@click.option("--medium_disk", "medium_disk", type=int, default=100)
@click.option("--huge_disk", "huge_disk", type=int, default=200)
@click.option(
    "--snp_recalibration_tranche_values",
    "snp_recalibration_tranche_values",
    multiple=True,
    type=str,
    default=[
        "100.0",
        "99.95",
        "99.9",
        "99.8",
        "99.6",
        "99.5",
        "99.4",
        "99.3",
        "99.0",
        "98.0",
        "97.0",
        "90.0",
    ],
)
@click.option(
    "--snp_recalibration_annotation_values",
    "snp_recalibration_annotation_values",
    multiple=True,
    type=str,
    default=[
        "AS_QD",
        "AS_MQRankSum",
        "AS_ReadPosRankSum",
        "AS_FS",
        "AS_SOR",
        "AS_MQ",
    ],
)
@click.option(
    "--indel_recalibration_tranche_values",
    "indel_recalibration_tranche_values",
    multiple=True,
    type=str,
    default=[
        "100.0",
        "99.95",
        "99.9",
        "99.5",
        "99.0",
        "97.0",
        "96.0",
        "95.0",
        "94.0",
        "93.5",
        "93.0",
        "92.0",
        "91.0",
        "90.0",
    ],
)
@click.option(
    "--indel_recalibration_annotation_values",
    "indel_recalibration_annotation_values",
    multiple=True,
    type=str,
    default=["AS_FS", "AS_SOR", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD"],
)
@click.option(
    "--eval_interval_list",
    "eval_interval_list",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list",
)
@click.option(
    "--hapmap_resource_vcf",
    "hapmap_resource_vcf",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz",
)
@click.option(
    "--hapmap_resource_vcf_index",
    "hapmap_resource_vcf_index",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi",
)
@click.option(
    "--omni_resource_vcf",
    "omni_resource_vcf",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz",
)
@click.option(
    "--omni_resource_vcf_index",
    "omni_resource_vcf_index",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi",
)
@click.option(
    "--one_thousand_genomes_resource_vcf",
    "one_thousand_genomes_resource_vcf",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
)
@click.option(
    "--one_thousand_genomes_resource_vcf_index",
    "one_thousand_genomes_resource_vcf_index",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
)
@click.option(
    "--mills_resource_vcf",
    "mills_resource_vcf",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
)
@click.option(
    "--mills_resource_vcf_index",
    "mills_resource_vcf_index",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
)
@click.option(
    "--axiom_poly_resource_vcf",
    "axiom_poly_resource_vcf",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
)
@click.option(
    "--axiom_poly_resource_vcf_index",
    "axiom_poly_resource_vcf_index",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi",
)
@click.option("--dbsnp_resource_vcf", "dbsnp_resource_vcf", type=str)
@click.option("--dbsnp_resource_vcf_index", "dbsnp_resource_vcf_index", type=str)
# ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
# than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
@click.option("--excess_het_threshold", "excess_het_threshold", type=float, default=54.69)
@click.option("--snp_filter_level", "snp_filter_level", type=float, default=99.7)
@click.option("--indel_filter_level", "indel_filter_level", type=float, default=99.0)
@click.option("--snp_vqsr_downsample_factor", "snp_vqsr_downsample_factor", type=int, default=10)
@click.option("--skip_allele_specific_annotations", "skip_allele_specific_annotations", is_flag=True)
@click.option("--dry_run", "dry_run", is_flag=True, default=False)
@click.option("--keep_scratch", "keep_scratch", is_flag=True, default=False)
@click.option("--is_small_callset", "is_small_callset", is_flag=True, default=False)
@click.option("--billing_project", "billing_project", type=str)
def main(
    combined_gvcf: str,
    output_bucket: str,
    num_gvcfs: int,
    callset_name: str,
    unpadded_intervals_file: str,
    ref_fasta: str,
    ref_fasta_index: str,
    ref_dict: str,
    dbsnp_vcf: str,
    dbsnp_vcf_index: str,
    small_disk: Optional[int] = None,
    medium_disk: Optional[int] = None,
    huge_disk: Optional[int] = None,
    snp_recalibration_tranche_values: List[str] = None,
    snp_recalibration_annotation_values: List[str] = None,
    indel_recalibration_tranche_values: List[str] = None,
    indel_recalibration_annotation_values: List[str] = None,
    eval_interval_list: str = None,
    hapmap_resource_vcf: str = None,
    hapmap_resource_vcf_index: str = None,
    omni_resource_vcf: str = None,
    omni_resource_vcf_index: str = None,
    one_thousand_genomes_resource_vcf: str = None,
    one_thousand_genomes_resource_vcf_index: str = None,
    mills_resource_vcf: str = None,
    mills_resource_vcf_index: str = None,
    axiom_poly_resource_vcf: str = None,
    axiom_poly_resource_vcf_index: str = None,
    dbsnp_resource_vcf: str = None,
    dbsnp_resource_vcf_index: str = None,
    excess_het_threshold: float = None,
    snp_filter_level: float = None,
    indel_filter_level: float = None,
    snp_vqsr_downsample_factor: int = None,
    skip_allele_specific_annotations: bool = None,
    dry_run: bool = None,
    keep_scratch: bool = None,
    billing_project: str = None,
    is_small_callset: bool = None,
):
    if not dry_run:
        if not billing_project:
            raise click.BadParameter(
                "--billing_project has to be specified (unless --dry_run is set)"
            )

    # Make a 2.5:1 interval number to samples in callset ratio interval list.
    # We allow overriding the behavior by specifying the desired number of vcfs
    # to scatter over for testing / special requests.
    scatter_count_scale_factor = 0.15
    scatter_count = int(round(scatter_count_scale_factor * num_gvcfs))
    scatter_count = max(scatter_count, 2)

    # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
    # We allow overriding this default behavior for testing / special requests.
    is_small_callset = is_small_callset if is_small_callset is not None else num_gvcfs < 1000

    backend = hb.ServiceBackend(billing_project=billing_project)
    b = hb.Batch("VariantCallingOFTHEFUTURE", backend=backend)

    combined_gvcf = b.read_input(combined_gvcf)
    ref_fasta = b.read_input_group(
        base=ref_fasta,
        dict=ref_dict or (
            ref_fasta.replace(".fasta", "")
            .replace(".fna", "")
            .replace(".fa", "") + ".dict"),
        fai=ref_fasta_index or (ref_fasta + ".fai"),
    )
    dbsnp_vcf = b.read_input_group(
        base=dbsnp_vcf, 
        index=dbsnp_vcf_index)
    eval_interval_list = b.read_input(eval_interval_list)
    hapmap_resource_vcf = b.read_input_group(
        base=hapmap_resource_vcf, 
        index=hapmap_resource_vcf_index)
    omni_resource_vcf = b.read_input_group(
        base=omni_resource_vcf, 
        index=omni_resource_vcf_index)
    one_thousand_genomes_resource_vcf = b.read_input_group(
        base=one_thousand_genomes_resource_vcf, 
        index=one_thousand_genomes_resource_vcf_index)
    mills_resource_vcf = b.read_input_group(
        base=mills_resource_vcf,
        index=mills_resource_vcf_index)
    axiom_poly_resource_vcf = b.read_input_group(
        base=axiom_poly_resource_vcf, 
        index=axiom_poly_resource_vcf_index)
    dbsnp_resource_vcf = b.read_input_group(
        base=dbsnp_resource_vcf,
        index=dbsnp_resource_vcf_index) if dbsnp_resource_vcf else dbsnp_vcf

    intervals = add_SplitIntervals_step(
        b,
        unpadded_intervals_file,
        scatter_count,
        ref_fasta,
        disk_size=small_disk,
    ).intervals

    gnarly_output_vcfs = [
        add_GnarlyGenotyperOnVcf_step(
            b,
            combined_gvcf=combined_gvcf,
            interval=intervals[f'interval_{idx}'],
            ref_fasta=ref_fasta,
            dbsnp_vcf=dbsnp_vcf,
            disk_size=medium_disk,
        ).output_vcf for idx in range(scatter_count)]

    if not is_small_callset:
        # ExcessHet filtering applies only to callsets with a large number of samples, 
        # e.g. hundreds of unrelated samples. Small cohorts should not trigger ExcessHet 
        # filtering as values should remain small. Note cohorts of consanguinous samples 
        # will inflate ExcessHet, and it is possible to limit the annotation to founders 
        # for such cohorts by providing a pedigree file during variant calling.
        hard_filtered_vcfs = [
            add_HardFilter_step(
                b,
                input_vcf=gnarly_output_vcfs[idx],
                excess_het_threshold=excess_het_threshold,
                disk_size=medium_disk,
            ).output_vcf
            for idx in range(scatter_count)]
    else:
        hard_filtered_vcfs = gnarly_output_vcfs
    # hard_filtered_vcfs = [
    #     b.read_input_group(
    #         base=f"gs://playground-au/batch/859e9a/{idx + 2}/output_vcf.vcf.gz", 
    #         index=f"gs://playground-au/batch/859e9a/{idx + 2}/output_vcf.vcf.gz.tbi"
    #     )
    #     for idx in range(scatter_count)
    # ]

    sites_only_vcfs = [
        add_MakeSitesOnlyVcf_step(
            b,
            input_vcf=hard_filtered_vcfs[idx],
            disk_size=medium_disk,
        ).sites_only_vcf
        for idx in range(scatter_count)]
    # sites_only_vcfs = [
    #     b.read_input_group(
    #         base=f"gs://playground-au/batch/859e9a/{idx + 9}/sites_only_vcf.vcf.gz", 
    #         index=f"gs://playground-au/batch/859e9a/{idx + 9}/sites_only_vcf.vcf.gz.tbi"
    #     )
    #     for idx in range(scatter_count)
    # ]

    sites_only_gathered_vcf = add_SitesOnlyGatherVcf_step(
        b,
        input_vcfs=sites_only_vcfs,
        disk_size=medium_disk,
    ).output_vcf
    
    indels_variant_recalibrator_job = add_IndelsVariantRecalibrator_step(
        b,
        sites_only_variant_filtered_vcf=sites_only_gathered_vcf,
        recalibration_tranche_values=indel_recalibration_tranche_values,
        recalibration_annotation_values=indel_recalibration_annotation_values,
        mills_resource_vcf=mills_resource_vcf,
        axiom_poly_resource_vcf=axiom_poly_resource_vcf,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        use_allele_specific_annotations=not skip_allele_specific_annotations,
        disk_size=small_disk,
    )
    indels_recalibration = indels_variant_recalibrator_job.recalibration
    indels_tranches = indels_variant_recalibrator_job.tranches
    # indels_recalibration = "gs://playground-au/batch/859e9a/17/recalibration"
    # indels_tranches = "gs://playground-au/batch/859e9a/17/tranches"

    model_report = add_SNPsVariantRecalibratorCreateModel_step(
        b,
        sites_only_variant_filtered_vcf=sites_only_gathered_vcf,
        recalibration_tranche_values=snp_recalibration_tranche_values,
        recalibration_annotation_values=snp_recalibration_annotation_values,
        downsample_factor=snp_vqsr_downsample_factor,
        hapmap_resource_vcf=hapmap_resource_vcf,
        omni_resource_vcf=omni_resource_vcf,
        one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        disk_size=small_disk,
        use_allele_specific_annotations=not skip_allele_specific_annotations,
        is_small_callset=is_small_callset,
    ).model_report
    # model_report = b.read_input("gs://playground-au/batch/859e9a/18/model_report")

    snps_recalibrator_jobs = [add_SNPsVariantRecalibratorScattered_step(
        b,
        sites_only_variant_filtered_vcf=sites_only_vcfs[idx],
        recalibration_tranche_values=snp_recalibration_tranche_values,
        recalibration_annotation_values=snp_recalibration_annotation_values,
        model_report=model_report,
        hapmap_resource_vcf=hapmap_resource_vcf,
        omni_resource_vcf=omni_resource_vcf,
        one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        disk_size=small_disk,
        is_small_callset=is_small_callset,
        use_allele_specific_annotations=not skip_allele_specific_annotations,
    ) for idx in range(len(sites_only_vcfs))]
    snps_recalibrations = [j.recalibrations for j in snps_recalibrator_jobs]
    snps_tranches = [j.tranches for j in snps_recalibrator_jobs]
    # snp_tranches = [
    #     b.read_input(f'gs://playground-au/batch/df311d/{idx + 1}/tranches') 
    #     for idx in range(scatter_count)
    # ]
    # snp_recalibrations = [
    #     b.read_input(f'gs://playground-au/batch/df311d/{idx + 1}/recalibration') 
    #     for idx in range(scatter_count)
    # ]

    snps_gathered_tranches = add_SNPGatherTranches_step(
        b,
        tranches=snps_tranches,
        disk_size=small_disk,
    ).out_tranches

    recalibrated_vcfs = [
        add_ApplyRecalibration_step(
            b,
            input_vcf=hard_filtered_vcfs[idx],
            indels_recalibration=indels_recalibration,
            indels_tranches=indels_tranches,
            snps_recalibration=snps_recalibrations[idx],
            snps_tranches=snps_tranches,
            indel_filter_level=indel_filter_level,
            snp_filter_level=snp_filter_level,
            disk_size=medium_disk,
            use_allele_specific_annotations=not skip_allele_specific_annotations,
        ).recalibrated_vcf for idx in range(len(hard_filtered_vcfs))
    ]
    
    final_gathered_vcf = add_FinalGatherVcf_step(
        b,
        input_vcfs=recalibrated_vcfs,
        output_vcf_path=os.path.join(output_bucket, callset_name + "-recalibrated.vcf.gz"),
        disk_size=huge_disk,
    ).output_vcf

    add_VariantEval_step(
        b,
        input_vcf=final_gathered_vcf,
        ref_fasta=ref_fasta,
        dbsnp_vcf=dbsnp_vcf,
        output_path=os.path.join(output_bucket, callset_name + "-cohorteval.txt"),
        disk_size=huge_disk,
    )

    b.run(dry_run=dry_run, delete_scratch_on_exit=not keep_scratch)


def add_SplitIntervals_step(
    b,
    interval_list,
    scatter_count,
    ref_fasta,
    disk_size,
):
    j = b.new_job("SplitIntervals")
    j.image(GATK_DOCKER)
    j.memory(f"8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(
        intervals={
            f"interval_{idx}":
            f"{{root}}/{str(idx).zfill(4)}-scattered.interval_list"
            for idx in range(scatter_count)
        })

    j.command(
        f"""set -e

    gatk --java-options -Xms3g SplitIntervals \\
      -L {interval_list} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {ref_fasta.base} \\
      -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
      """
    )
    return j


def add_GnarlyGenotyperOnVcf_step(
    b,
    combined_gvcf,
    interval,
    ref_fasta,
    dbsnp_vcf,
    disk_size,
):
    j = b.new_job("GnarlyGenotyperOnVcf")
    j.image(GNARLY_DOCKER)  # GnarlyGenotyper crashes with NullPointerException when using GATK docker
    j.memory(f"32G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(output_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"})

    j.command(
        f"""set -e

    tabix -p vcf {combined_gvcf}

    gatk --java-options -Xms8g \\
      GnarlyGenotyper \\
      -R {ref_fasta.base} \\
      -O {j.output_vcf["vcf.gz"]} \\
      -D {dbsnp_vcf.base} \\
      --only-output-calls-starting-in-intervals \\
      --keep-all-sites \\
      -V {combined_gvcf} \\
      -L {interval} \\
      --create-output-variant-index"""
    )
    return j


def add_HardFilter_step(
    b,
    input_vcf,
    excess_het_threshold,
    disk_size,
):
    j = b.new_job("HardFilter")
    j.image(GATK_DOCKER)
    j.memory("8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(output_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"})

    j.command(
        f"""set -euo pipefail

    # Captring stderr to avoid Batch pod from crashing with OOM from millions of warning messages
    # from VariantFiltration (eg. `JexlEngine - ![0,9]: 'ExcessHet > 54.69;' undefined variable ExcessHet`)
    gatk --java-options -Xms3g \\
      VariantFiltration \\
      --filter-expression 'ExcessHet > {excess_het_threshold}' \\
      --filter-name ExcessHet \\
      -O {j.output_vcf["vcf.gz"]} \\
      -V {input_vcf["vcf.gz"]} \\
      2> {j.stderr}
    """
    )
    return j


def add_MakeSitesOnlyVcf_step(
    b,
    input_vcf,
    disk_size,
):
    j = b.new_job("MakeSitesOnlyVcf")
    j.image(GATK_DOCKER)
    j.memory("8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(sites_only_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"})

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      MakeSitesOnlyVcf \\
      -I {input_vcf["vcf.gz"]} \\
      -O {j.sites_only_vcf["vcf.gz"]}"""
    )
    return j


def add_SitesOnlyGatherVcf_step(
    b,
    input_vcfs,
    disk_size,
):
    j = b.new_job("SitesOnlyGatherVcf")
    j.image(GATK_DOCKER)
    j.memory("8G")
    j.storage(f"{disk_size}G")

    j.declare_resource_group(output_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"})

    input_cmdl = " ".join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options -Xms6g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      {input_cmdl} \\
      --output {j.output_vcf["vcf.gz"]}

    tabix {j.output_vcf["vcf.gz"]}"""
    )
    return j


def add_IndelsVariantRecalibrator_step(
    b,
    recalibration_tranche_values,
    recalibration_annotation_values,
    sites_only_variant_filtered_vcf,
    mills_resource_vcf,
    axiom_poly_resource_vcf,
    dbsnp_resource_vcf,
    use_allele_specific_annotations,
    disk_size,
    max_gaussians=4,
):
    j = b.new_job("IndelsVariantRecalibrator")
    j.image(GATK_DOCKER)
    j.memory(f"32G")
    j.cpu(2)
    j.storage(f"{disk_size}G")

    j.declare_resource_group(recalibration={"index": "{root}.idx", "base": "{root}"})
    
    tranche_cmdl = " ".join([f"-tranche {v}" for v in recalibration_tranche_values])
    an_cmdl = " ".join([f"-an {v}" for v in recalibration_annotation_values])
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms24g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf["vcf.gz"]} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode INDEL \\
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''} \\
      --max-gaussians {max_gaussians} \\
      -resource:mills,known=false,training=true,truth=true,prior=12 {mills_resource_vcf.base} \\
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {axiom_poly_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 {dbsnp_resource_vcf.base} \\
      --rscript-file {j.indel_rscript_file}"""
    )
    return j


def add_SNPsVariantRecalibratorCreateModel_step(
    b,
    downsample_factor,
    recalibration_tranche_values,
    recalibration_annotation_values,
    sites_only_variant_filtered_vcf,
    hapmap_resource_vcf,
    omni_resource_vcf,
    one_thousand_genomes_resource_vcf,
    dbsnp_resource_vcf,
    disk_size,
    use_allele_specific_annotations,
    is_small_callset=False,
    max_gaussians=None,
):
    j = b.new_job("SNPsVariantRecalibratorCreateModel")
    j.image(GATK_DOCKER)
    mem_gb = 64 if is_small_callset else 128
    j.memory(f"{mem_gb}G")
    j.cpu(2)
    j.storage(f"{disk_size}G")

    max_gaussians = 4 if is_small_callset else 6

    tranche_cmdl = " ".join([f"-tranche {v}" for v in recalibration_tranche_values])
    an_cmdl = " ".join([f"-an {v}" for v in recalibration_annotation_values])
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{mem_gb - 2}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf["vcf.gz"]} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''} \\
      --sample-every-Nth-variant {downsample_factor} \\
      --output-model {j.model_report} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base} \\
      --rscript-file {j.snp_rscript_file}"""
    )
    return j


def add_SNPsVariantRecalibratorScattered_step(
    b,
    sites_only_variant_filtered_vcf,
    recalibration_tranche_values,
    recalibration_annotation_values,
    hapmap_resource_vcf,
    omni_resource_vcf,
    one_thousand_genomes_resource_vcf,
    dbsnp_resource_vcf,
    disk_size,
    use_allele_specific_annotations,
    is_small_callset=False,
    model_report=None,
    max_gaussians=None,
):
    j = b.new_job("SNPsVariantRecalibratorScattered")
    
    mem_gb = 32 if is_small_callset else 64

    model_report_arg = (
        f"--input-model {model_report} --output-tranches-for-scatter"
        if model_report is not None
        else ""
    )
    j.image(GATK_DOCKER)
    j.memory(f"{mem_gb}G")
    j.cpu(2)
    j.storage(f"{disk_size}G")

    j.declare_resource_group(recalibration={"index": "{root}.idx", "base": "{root}"})
    
    max_gaussians = 4 if is_small_callset else 6

    tranche_cmdl = " ".join([f"-tranche {v}" for v in recalibration_tranche_values])
    an_cmdl = " ".join([f"-an {v}" for v in recalibration_annotation_values])
    j.command(
        f"""set -euo pipefail

    MODEL_REPORT={model_report}

    gatk --java-options -Xms{mem_gb - 1}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf["vcf.gz"]} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''} \\
      {model_report_arg} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base}"""
    )
    return j


def add_SNPGatherTranches_step(
    b,
    tranches,
    disk_size,
):
    j = b.new_job("SNPGatherTranches")
    j.image(GATK_DOCKER)
    j.memory("8G")
    j.cpu(2)
    j.storage(f"{disk_size}G")

    inputs_cmdl = " ".join([f"--input {t}" for t in tranches])
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      GatherTranches \\
      --mode SNP \\
      {inputs_cmdl} \\
      --output {j.out_tranches}"""
    )
    return j


def add_ApplyRecalibration_step(
    b,
    input_vcf,
    indels_recalibration,
    indels_tranches,
    snps_recalibration,
    snps_tranches,
    indel_filter_level,
    snp_filter_level,
    use_allele_specific_annotations,
    disk_size,
):
    j = b.new_job("ApplyRecalibration")
    j.image(GATK_DOCKER)
    j.memory("8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(recalibrated_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"})

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O tmp.indel.recalibrated.vcf \\
      -V {input_vcf["vcf.gz"]} \\
      --recal-file {indels_recalibration} \\
      --tranches-file {indels_tranches} \\
      --truth-sensitivity-filter-level {indel_filter_level} \\
      --create-output-variant-index true \\
      -mode INDEL \\
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''}

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O {j.recalibrated_vcf["vcf.gz"]} \\
      -V tmp.indel.recalibrated.vcf \\
      --recal-file {snps_recalibration} \\
      --tranches-file {snps_tranches} \\
      --truth-sensitivity-filter-level {snp_filter_level} \\
      --create-output-variant-index true \\
      -mode SNP \\
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''}"""
    )
    return j


def add_CollectMetricsSharded_step(
    b,
    input_vcf,
    dbsnp_vcf,
    interval_list,
    ref_dict,
    disk_size,
):
    j = b.new_job("CollectMetricsSharded")
    j.image(GATK_DOCKER)
    j.memory("8G")
    j.cpu(2)
    j.storage(f"{disk_size}G")
    j.declare_resource_group(metrics={
        "detail_metrics": "{root}.variant_calling_detail_metrics",
        "summary_metrics": "{root}.variant_calling_summary_metrics",
    })

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      CollectVariantCallingMetrics \\
      --INPUT {input_vcf["vcf.gz"]} \\
      --DBSNP {dbsnp_vcf.base} \\
      --SEQUENCE_DICTIONARY {ref_dict} \\
      --OUTPUT {j.metrics} \\
      --THREAD_COUNT 8 \\
      --TARGET_INTERVALS {interval_list}"""
    )
    return j


def add_FinalGatherVcf_step(
    b,
    input_vcfs,
    output_vcf_path,
    disk_size,
):
    j = b.new_job("FinalGatherVcf")
    j.image(GATK_DOCKER)
    j.memory(f"8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(output_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"})

    input_cmdl = " ".join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options -Xms6g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      {input_cmdl} \\
      --output {j.output_vcf["vcf.gz"]}

    tabix {j.output_vcf}"""
    )
    b.write_output(j.output_vcf, output_vcf_path.replace(".vcf.gz", ""))
    return j


def add_VariantEval_step(
    b,
    input_vcf,
    ref_fasta,
    dbsnp_vcf,
    output_path,
    disk_size,
):
    j = b.new_job("VariantEval")
    j.image(GATK_DOCKER)
    j.memory(f"8G")
    j.storage(f"{disk_size}G")

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      VariantEval \\
      --eval {input_vcf} \\
      -R {ref_fasta.base} \\
      -D {dbsnp_vcf.base} \\
      --output {j.output}"""
    )
    b.write_output(j.output, output_path)
    return j


def add_GatherVariantCallingMetrics_step(
    b,
    input_details,
    input_summaries,
    output_path_prefix,
    disk_size,
):
    j = b.new_job("GatherVariantCallingMetrics")
    j.image(GATK_DOCKER)
    j.memory(f"8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(metrics={
        "detail_metrics": "{root}.variant_calling_detail_metrics",
        "summary_metrics": "{root}.variant_calling_summary_metrics",
    })

    input_cmdl = " ".join("--INPUT {f} " for f in input_details + input_summaries)
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms2g \\
      AccumulateVariantCallingMetrics \\
      {input_cmdl} \\
      --OUTPUT {j.metrics}"""
    )
    b.write_output(j.metrics, output_path_prefix)
    return j


if __name__ == "__main__":
    main()
