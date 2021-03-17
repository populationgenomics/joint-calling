import math, os
from typing import Union, Optional, List

import click

import hailtop.batch as hb


@click.command()
@click.option("--combined_gvcf", "combined_gvcf", type=str, required=True)
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
@click.option("--small_disk", "small_disk", type=int, default=100)
@click.option("--medium_disk", "medium_disk", type=int, default=200)
@click.option("--huge_disk", "huge_disk", type=int, default=2000)
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
@click.option("--indel_vqsr_downsample_factor", "indel_vqsr_downsample_factor", type=int, default=10)
@click.option("--skip_allele_specific_annotations", "skip_allele_specific_annotations", is_flag=True)
@click.option("--dry_run", "dry_run", is_flag=True, default=False)
@click.option("--is_small_callset", "is_small_callset", is_flag=True, default=False)
@click.option("--billing_project", "billing_project", type=str)
def main(
    combined_gvcf: str,
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
    indel_vqsr_downsample_factor: int = None,
    skip_allele_specific_annotations: bool = None,
    dry_run: bool = None,
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
    # Zamboni notes say "WGS runs get 30x more scattering than Exome" and
    # exome scatterCountPerSample is 0.05, min scatter 10, max 1000
    # 
    # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
    # We allow overriding this default behavior for testing / special requests.
    is_small_callset = is_small_callset if is_small_callset is not None else num_gvcfs < 1000
    scatter_count_scale_factor = 0.15
    scatter_count = int(round(scatter_count_scale_factor * num_gvcfs))
    scatter_count = max(scatter_count, 2)

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
        tbi=dbsnp_vcf_index)
    eval_interval_list = b.read_input(eval_interval_list)
    hapmap_resource_vcf = b.read_input_group(
        base=hapmap_resource_vcf, 
        tbi=hapmap_resource_vcf_index)
    omni_resource_vcf = b.read_input_group(
        base=omni_resource_vcf, 
        tbi=omni_resource_vcf_index)
    one_thousand_genomes_resource_vcf = b.read_input_group(
        base=one_thousand_genomes_resource_vcf, 
        tbi=one_thousand_genomes_resource_vcf_index)
    mills_resource_vcf = b.read_input_group(
        base=mills_resource_vcf,
        tbi=mills_resource_vcf_index)
    axiom_poly_resource_vcf = b.read_input_group(
        base=axiom_poly_resource_vcf, 
        tbi=axiom_poly_resource_vcf_index)
    dbsnp_resource_vcf = b.read_input_group(
        base=dbsnp_resource_vcf,
        tbi=dbsnp_resource_vcf_index) if dbsnp_resource_vcf else dbsnp_vcf

    SplitIntervals = add_SplitIntervals_step(
        b,
        unpadded_intervals_file,
        scatter_count,
        ref_fasta,
        disk_size=small_disk,
    )

    GnarlyGenotyperOnVcf = []
    for idx in range(scatter_count):
        GnarlyGenotyperOnVcf.append(
            add_GnarlyGenotyperOnVcf_step(
                b,
                combined_gvcf=combined_gvcf,
                interval=SplitIntervals.intervals[f'interval_{idx}'],
                output_vcf_filename=f"{callset_name}.{idx}.vcf.gz",
                ref_fasta=ref_fasta,
                dbsnp_vcf=dbsnp_vcf,
                is_small_callset=is_small_callset,
            )
        )

    HardFilterAndMakeSitesOnlyVcf = []
    for idx in range(scatter_count):
        HardFilterAndMakeSitesOnlyVcf.append(
            add_HardFilterAndMakeSitesOnlyVcf_step(
                b,
                input_vcf=GnarlyGenotyperOnVcf[idx].output_vcf,
                excess_het_threshold=excess_het_threshold,
                variant_filtered_vcf_filename=f"{callset_name}.{idx}.variant_filtered.vcf.gz",
                sites_only_vcf_filename=f"{callset_name}.{idx}.sites_only.variant_filtered.vcf.gz",
                disk_size=medium_disk,
            )
        )
    SitesOnlyGatherVcf = add_SitesOnlyGatherVcf_step(
        b,
        input_vcfs=[j.sites_only_vcf for j in HardFilterAndMakeSitesOnlyVcf],
        output_vcf_name=f"{callset_name}.sites_only.vcf.gz",
        disk_size=medium_disk,
    )
    IndelsVariantRecalibrator = add_IndelsVariantRecalibrator_step(
        b,
        sites_only_variant_filtered_vcf=SitesOnlyGatherVcf.output_vcf,
        recalibration_filename=f"{callset_name}.indels.recal",
        tranches_filename=f"{callset_name}.indels.tranches",
        recalibration_tranche_values=indel_recalibration_tranche_values,
        recalibration_annotation_values=indel_recalibration_annotation_values,
        mills_resource_vcf=mills_resource_vcf,
        axiom_poly_resource_vcf=axiom_poly_resource_vcf,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        use_allele_specific_annotations=not skip_allele_specific_annotations,
        disk_size=small_disk,
    )
    SNPsVariantRecalibratorCreateModel = add_SNPsVariantRecalibratorCreateModel_step(
        b,
        sites_only_variant_filtered_vcf=SitesOnlyGatherVcf.output_vcf,
        recalibration_filename=f"{callset_name}.snps.recal",
        tranches_filename=f"{callset_name}.snps.tranches",
        recalibration_tranche_values=snp_recalibration_tranche_values,
        recalibration_annotation_values=snp_recalibration_annotation_values,
        downsample_factor=snp_vqsr_downsample_factor,
        model_report_filename=f"{callset_name}.snps.model.report",
        hapmap_resource_vcf=hapmap_resource_vcf,
        omni_resource_vcf=omni_resource_vcf,
        one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        disk_size=small_disk,
        use_allele_specific_annotations=not skip_allele_specific_annotations,
        is_small_callset=is_small_callset,
    )

    SNPsVariantRecalibratorScattered = []
    for j in HardFilterAndMakeSitesOnlyVcf:
        SNPsVariantRecalibratorScattered.append(
            add_SNPsVariantRecalibratorScattered_step(
                b,
                sites_only_variant_filtered_vcf=j.sites_only_vcf,
                recalibration_filename=f"{callset_name}.snps.{idx}.recal",
                tranches_filename=f"{callset_name}.snps.{idx}.tranches",
                recalibration_tranche_values=snp_recalibration_tranche_values,
                recalibration_annotation_values=snp_recalibration_annotation_values,
                model_report=SNPsVariantRecalibratorCreateModel.model_report,
                hapmap_resource_vcf=hapmap_resource_vcf,
                omni_resource_vcf=omni_resource_vcf,
                one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
                dbsnp_resource_vcf=dbsnp_resource_vcf,
                disk_size=small_disk,
                is_small_callset=is_small_callset,
                use_allele_specific_annotations=not skip_allele_specific_annotations,
            )
        )
    SNPGatherTranches = add_SNPGatherTranches_step(
        b,
        tranches=[j.tranches for j in SNPsVariantRecalibratorScattered],
        output_filename=f"{callset_name}.snps.gathered.tranches",
        disk_size=small_disk,
    )
    # Need to call `depends_on` explicitly because doesn't resolve the dependency from the tranches input
    SNPGatherTranches.depends_on(*SNPsVariantRecalibratorScattered)

    ApplyRecalibration = []
    for idx, j in enumerate(HardFilterAndMakeSitesOnlyVcf):
        ApplyRecalibration.append(
            add_ApplyRecalibration_step(
                b,
                recalibrated_vcf_filename=f"{callset_name}.filtered.{idx}.vcf.gz",
                input_vcf=j.variant_filtered_vcf,
                indels_recalibration=IndelsVariantRecalibrator.recalibration,
                indels_tranches=IndelsVariantRecalibrator.tranches,
                snps_recalibration=SNPsVariantRecalibratorScattered[idx].recalibration,
                snps_tranches=SNPGatherTranches.out_tranches,
                indel_filter_level=indel_filter_level,
                snp_filter_level=snp_filter_level,
                disk_size=medium_disk,
                use_allele_specific_annotations=not skip_allele_specific_annotations,
            )
        )
    # For large callsets we need to collect metrics from the shards and gather them later
    CollectMetricsSharded = None
    if not is_small_callset:
        CollectMetricsSharded = add_CollectMetricsSharded_step(
            b,
            input_vcf=ApplyRecalibration.recalibrated_vcf,
            metrics_filename_prefix=f"{callset_name}.{idx}",
            dbsnp_vcf=dbsnp_vcf,
            interval_list=eval_interval_list,
            ref_dict=ref_fasta.dict,
            disk_size=small_disk,
        )

    # For small callsets we can gather the VCF shards and then collect metrics on it
    FinalGatherVcf = None
    CollectMetricsOnFullVcf = None
    if is_small_callset:
        FinalGatherVcf = add_FinalGatherVcf_step(
            b,
            input_vcfs=[j.recalibrated_vcf for j in ApplyRecalibration],
            output_vcf_name=f"{callset_name}.vcf.gz",
            disk_size=huge_disk,
        )
    
        CollectMetricsOnFullVcf = add_CollectMetricsOnFullVcf_step(
            b,
            input_vcf=FinalGatherVcf.output_vcf,
            metrics_filename_prefix=callset_name,
            dbsnp_vcf=dbsnp_vcf,
            interval_list=eval_interval_list,
            ref_dict=ref_fasta.dict,
            disk_size=huge_disk,
        )
  
    # For large callsets we still need to gather the sharded metrics
    if not is_small_callset:
        GatherVariantCallingMetrics = add_GatherVariantCallingMetrics_step(
            b,
            input_details = [j.detail_metrics_file for j in CollectMetricsSharded],
            input_summaries = [j.summary_metrics_file for j in CollectMetricsSharded],
            output_prefix = callset_name,
            disk_size = medium_disk
        )
    
    b.run(dry_run=dry_run)


def add_GnarlyGenotyperOnVcf_step(
    b,
    combined_gvcf,
    interval,
    output_vcf_filename,
    ref_fasta,
    dbsnp_vcf,
    gatk_docker="gcr.io/broad-dsde-methods/gnarly_genotyper:hail_ukbb_300K",
    is_small_callset=False,
    disk_size=None,
):
    j = b.new_job("GnarlyGenotyperOnVcf")
    disk_size = (
        disk_size if disk_size is not None else (40 if is_small_callset else 80)
    )
    j.declare_resource_group(output_vcf={"index": "{root}.vcf.gz.tbi", "base": "{root}.vcf.gz"})
    j.image(gatk_docker)
    j.memory(f"32G")
    j.storage(f"{disk_size}G")

    j.command(
        f"""set -e

    tabix -p vcf {combined_gvcf}

    gatk --java-options -Xms8g \\
      GnarlyGenotyper \\
      -R {ref_fasta.base} \\
      -O {j.output_vcf.base} \\
      -D {dbsnp_vcf.base} \\
      --only-output-calls-starting-in-intervals \\
      --keep-all-sites \\
      -V {combined_gvcf} \\
      -L {interval} \\
      --create-output-variant-index"""
    )
    b.write_output(j.output_vcf, output_vcf_filename)
    return j


def add_SplitIntervals_step(
    b,
    interval_list,
    scatter_count,
    ref_fasta,
    disk_size,
):
    j = b.new_job("SplitIntervals")
    j.declare_resource_group(
        intervals={
            f"interval_{idx}":
            f"{{root}}/{str(idx).zfill(4)}-scattered.interval_list"
            for idx in range(scatter_count)
        })
    j.image("us.gcr.io/broad-gatk/gatk:4.1.8.0")
    j.memory(f"8G")
    j.storage(f"{disk_size}G")

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
    b.write_output(j.intervals, "scatter_intervals")
    # for scatter_idx in range(scatter_count):
    #     interval_file_prefix = str(scatter_idx).zfill(4)
    #     b.write_output(j.scatter_dir, f"{interval_file_prefix}-scattered.interval_list")
    return j


def add_HardFilterAndMakeSitesOnlyVcf_step(
    b,
    input_vcf,
    excess_het_threshold,
    variant_filtered_vcf_filename,
    sites_only_vcf_filename,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.8.0",
):
    j = b.new_job("HardFilterAndMakeSitesOnlyVcf")
    j.declare_resource_group(variant_filtered_vcf={"index": "{root}.vcf.gz.tbi", "base": "{root}.vcf.gz"})
    j.declare_resource_group(sites_only_vcf={"index": "{root}.vcf.gz.tbi", "base": "{root}.vcf.gz"})
    j.image(gatk_docker)
    j.memory("8G")
    j.storage(f"{disk_size}G")

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms3g \\
      VariantFiltration \\
      --filter-expression 'ExcessHet > {excess_het_threshold}' \\
      --filter-name ExcessHet \\
      -O {j.variant_filtered_vcf.base} \\
      -V {input_vcf.base}

    gatk --java-options -Xms3g \\
      MakeSitesOnlyVcf \\
      -I {j.variant_filtered_vcf.base} \\
      -O {j.sites_only_vcf.base}"""
    )
    b.write_output(j.variant_filtered_vcf, variant_filtered_vcf_filename)
    b.write_output(j.sites_only_vcf, sites_only_vcf_filename)
    return j


def add_SitesOnlyGatherVcf_step(
    b,
    input_vcfs,
    output_vcf_name,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.8.0",
):
    j = b.new_job("SitesOnlyGatherVcf")
    j.image(gatk_docker)
    j.memory("8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(output_vcf={"index": "{root}.vcf.gz.tbi", "base": "{root}.vcf.gz"})

    input_cmdl = " ".join([f"--input {v.base}" for v in input_vcfs])
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
      --output {j.output_vcf.base}

    tabix {j.output_vcf.base}"""
    )
    b.write_output(j.output_vcf, output_vcf_name)
    return j


def add_IndelsVariantRecalibrator_step(
    b,
    recalibration_filename,
    tranches_filename,
    recalibration_tranche_values,
    recalibration_annotation_values,
    sites_only_variant_filtered_vcf,
    mills_resource_vcf,
    axiom_poly_resource_vcf,
    dbsnp_resource_vcf,
    use_allele_specific_annotations,
    disk_size,
    max_gaussians=4,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.8.0",
):
    j = b.new_job("IndelsVariantRecalibrator")
    j.image(gatk_docker)
    j.memory(f"32G")
    j.cpu(2)
    j.storage(f"{disk_size}G")
    
    tranche_cmdl = " ".join([f"-tranche {v}" for v in recalibration_tranche_values])
    an_cmdl = " ".join([f"-an {v}" for v in recalibration_annotation_values])
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms24g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf.base} \\
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
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 {dbsnp_resource_vcf.base}"""
    )
    b.write_output(j.recalibration, recalibration_filename)
    b.write_output(j.tranches, tranches_filename)
    return j


def add_SNPsVariantRecalibratorCreateModel_step(
    b,
    recalibration_filename,
    tranches_filename,
    downsample_factor,
    model_report_filename,
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
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.4.1",
):
    j = b.new_job("SNPsVariantRecalibratorCreateModel")
    j.image(gatk_docker)
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
      -V {sites_only_variant_filtered_vcf.base} \\
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
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base}"""
    )
    b.write_output(j.recalibration, recalibration_filename)
    b.write_output(j.tranches, tranches_filename)
    b.write_output(j.model_report, model_report_filename)
    return j


def add_SNPsVariantRecalibratorScattered_step(
    b,
    recalibration_filename,
    tranches_filename,
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
    model_report=None,
    max_gaussians=None,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.8.0",
    model_report_arg=None,
):
    j = b.new_job("SNPsVariantRecalibratorScattered")
    
    mem_gb = 32 if is_small_callset else 64

    model_report_arg = (
        model_report_arg
        if model_report_arg is not None
        else (
            (("--input-model " + model_report) + "--output-tranches-for-scatter")
            if model_report is not None
            else ""
        )
    )
    j.image(gatk_docker)
    j.memory(f"{mem_gb}G")
    j.cpu(2)
    j.storage(f"{disk_size}G")

    max_gaussians = 4 if is_small_callset else 6

    tranche_cmdl = " ".join([f"-tranche {v}" for v in recalibration_tranche_values])
    j.command(
        f"""set -euo pipefail

    MODEL_REPORT={model_report}

    gatk --java-options -Xms{mem_gb - 1}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf.base} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      -an {recalibration_annotation_values} \\
      -mode SNP \\
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''} \\
      {model_report_arg} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base}"""
    )
    b.write_output(j.recalibration, recalibration_filename)
    b.write_output(j.tranches, tranches_filename)
    return j


def add_SNPGatherTranches_step(
    b,
    tranches,
    output_filename,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.8.0",
):
    j = b.new_job("SNPGatherTranches")
    j.image(gatk_docker)
    j.memory("8G")
    j.cpu(2)
    j.storage(f"{disk_size}G")

    tranches_file = "tranches.list"
    with open(tranches_file, 'w') as out:
        for tranch in tranches:
            out.write(f'{tranch}\n')

    j.command(
        f"""set -euo pipefail

    tranches_fofn={tranches_file}

    # Jose says:
    # Cromwell will fall over if we have it try to localize tens of thousands of files,
    # so we manually localize files using gsutil.
    # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
    # PAPI doesn't do.

    # This is here to deal with the JES bug where commands may be run twice
    rm -rf tranches
    mkdir tranches
    RETRY_LIMIT=5

    count=0
    until cat $tranches_fofn | /usr/bin/gsutil -m cp -L cp.log -c -I tranches/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ '$count' -ge '$RETRY_LIMIT' ]; then
        echo 'Could not copy all the tranches from the cloud' && exit 1
    fi

    cat $tranches_fofn | rev | cut -d '/' -f 1 | rev | awk '{{print 'tranches/' $1}}' > inputs.list

    gatk --java-options -Xms6g \\
      GatherTranches \\
      --input inputs.list \\
      --output {j.out_tranches}"""
    )
    b.write_output(j.out_tranches, output_filename)
    return j


def add_ApplyRecalibration_step(
    b,
    recalibrated_vcf_filename,
    input_vcf,
    indels_recalibration,
    indels_tranches,
    snps_recalibration,
    snps_tranches,
    indel_filter_level,
    snp_filter_level,
    use_allele_specific_annotations,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.8.0",
):
    j = b.new_job("ApplyRecalibration")
    j.image(gatk_docker)
    j.memory("8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(recalibrated_vcf={"index": "{root}.vcf.gz.tbi", "base": "{root}.vcf.gz"})

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O tmp.indel.recalibrated.vcf \\
      -V {input_vcf.base} \\
      --recal-file {indels_recalibration} \\
      --tranches-file {indels_tranches} \\
      --truth-sensitivity-filter-level {indel_filter_level} \\
      --create-output-variant-index true \\
      -mode INDEL \\
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''} \\

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O {j.recalibrated_vcf.base} \\
      -V tmp.indel.recalibrated.vcf \\
      --recal-file {snps_recalibration} \\
      --tranches-file {snps_tranches} \\
      --truth-sensitivity-filter-level {snp_filter_level} \\
      --create-output-variant-index true \\
      -mode SNP \\
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''} \\
"""
    )
    b.write_output(j.recalibrated_vcf, recalibrated_vcf_filename)
    return j


def add_CollectMetricsSharded_step(
    b,
    input_vcf,
    metrics_filename_prefix,
    dbsnp_vcf,
    interval_list,
    ref_dict,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.8.0",
):
    j = b.new_job("CollectMetricsSharded")
    j.image(gatk_docker)
    j.memory("8G")
    j.cpu(2)
    j.storage(f"{disk_size}G")
    j.declare_resource_group(metrics={
        "detail_metrics_file": "{root}.variant_calling_detail_metrics",
        "summary_metrics_file": "{root}.variant_calling_summary_metrics",
    })
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      CollectVariantCallingMetrics \\
      --INPUT {input_vcf.base} \\
      --DBSNP {dbsnp_vcf.base} \\
      --SEQUENCE_DICTIONARY {ref_dict} \\
      --OUTPUT {j.metrics} \\
      --THREAD_COUNT 8 \\
      --TARGET_INTERVALS {interval_list}"""
    )
    b.write_output(j.metrics, metrics_filename_prefix)
    return j


def add_FinalGatherVcf_step(
    b,
    input_vcfs,
    output_vcf_name,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.8.0",
):
    j = b.new_job("FinalGatherVcf")
    j.image(gatk_docker)
    j.memory(f"8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(output_vcf={"index": "{root}.vcf.gz.tbi", "base": "{root}.vcf.gz"})

    input_cmdl = " ".join([f"--input {v.base}" for v in input_vcfs])
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
      --output {j.output_vcf.base}

    tabix {j.output_vcf.base}"""
    )
    b.write_output(j.output_vcf, output_vcf_name)
    return j


def add_CollectMetricsOnFullVcf_step(
    b,
    input_vcf,
    metrics_filename_prefix,
    dbsnp_vcf,
    interval_list,
    ref_dict,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.8.0",
):
    j = b.new_job("CollectMetricsOnFullVcf")
    j.image(gatk_docker)
    j.memory(f"8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(metrics={
        "detail_metrics_file": "{root}.variant_calling_detail_metrics",
        "summary_metrics_file": "{root}.variant_calling_summary_metrics",
    })

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      CollectVariantCallingMetrics \\
      --INPUT {input_vcf.base} \\
      --DBSNP {dbsnp_vcf.base} \\
      --SEQUENCE_DICTIONARY {ref_dict} \\
      --OUTPUT {j.metrics} \\
      --THREAD_COUNT 8 \\
      --TARGET_INTERVALS {interval_list}"""
    )
    b.write_output(j.metrics, metrics_filename_prefix)
    return j


def add_GatherVariantCallingMetrics_step(
    b,
    input_details,
    input_summaries,
    output_prefix,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.8.0",
):
    j = b.new_job("GatherVariantCallingMetrics")
    j.image(gatk_docker)
    j.memory(f"8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(metrics={
        "detail_metrics_file": "{root}.variant_calling_detail_metrics",
        "summary_metrics_file": "{root}.variant_calling_summary_metrics",
    })

    input_cmdl = " ".join("--INPUT {f} " for f in input_details + input_summaries)

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms2g \\
      AccumulateVariantCallingMetrics \\
      {input_cmdl} \\
      --OUTPUT {j.metrics}"""
    )
    b.write_output(j.output, output_prefix)
    return j


def apply_secondary_file_format_to_filename(
    filepath: Optional[str], secondary_file: str
):
    """
    This is actually clever, you can probably trust this to do what you want.
    :param filepath: Filename to base
    :param secondary_file: CWL secondary format (Remove 1 extension for each leading ^.
    """
    if not filepath:
        return None

    fixed_sec = secondary_file.lstrip("^")
    leading = len(secondary_file) - len(fixed_sec)
    if leading <= 0:
        return filepath + fixed_sec

    basepath = ""
    filename = filepath
    if "/" in filename:
        idx = len(filepath) - filepath[::-1].index("/")
        basepath = filepath[:idx]
        filename = filepath[idx:]

    split = filename.split(".")

    newfname = filename + fixed_sec
    if len(split) > 1:
        newfname = ".".join(split[: -min(leading, len(split) - 1)]) + fixed_sec
    return basepath + newfname


if __name__ == "__main__":
    main()
