import math, os
from typing import Union, Optional, List

import click

import hailtop.batch as hb


def main(
    combined_gvcf: str,
    callset_name: str,
    unpadded_intervals_file: Optional[
        str
    ] = "gs://broad-references-private/HybSelOligos/xgen_plus_spikein/white_album_exome_calling_regions.v1.interval_list",
    ref_fasta: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
    ref_fasta_index: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
    ref_dict: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
    dbsnp_vcf: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
    dbsnp_vcf_index: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
    small_disk: Optional[int] = 100,
    medium_disk: Optional[int] = 200,
    huge_disk: Optional[int] = 2000,
    snp_recalibration_tranche_values: Optional[List[str]] = [
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
    snp_recalibration_annotation_values: Optional[List[str]] = [
        "AS_QD",
        "AS_MQRankSum",
        "AS_ReadPosRankSum",
        "AS_FS",
        "AS_SOR",
        "AS_MQ",
    ],
    indel_recalibration_tranche_values: Optional[List[str]] = [
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
    indel_recalibration_annotation_values: Optional[List[str]] = [
        "AS_FS",
        "AS_SOR",
        "AS_ReadPosRankSum",
        "AS_MQRankSum",
        "AS_QD",
    ],
    eval_interval_list: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/exome_evaluation_regions.v1.interval_list",
    hapmap_resource_vcf: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz",
    hapmap_resource_vcf_index: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi",
    omni_resource_vcf: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz",
    omni_resource_vcf_index: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi",
    one_thousand_genomes_resource_vcf: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
    one_thousand_genomes_resource_vcf_index: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
    mills_resource_vcf: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    mills_resource_vcf_index: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
    axiomPoly_resource_vcf: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
    axiomPoly_resource_vcf_index: Optional[
        str
    ] = "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi",
    dbsnp_resource_vcf: Optional[str] = None,
    dbsnp_resource_vcf_index: Optional[str] = None,
    excess_het_threshold: Optional[float] = 54.69,
    snp_filter_level: Optional[float] = 99.7,
    indel_filter_level: Optional[float] = 99.0,
    SNP_VQSR_downsampleFactor: Optional[int] = 75,
    indel_VQSR_downsampleFactor: Optional[int] = 10,
    use_allele_specific_annotations: Optional[bool] = True,
    vcf_count: Optional[int] = 30,
    unboundedScatterCount: Optional[int] = None,
    scatterCount: Optional[int] = None,
    unpadded_intervals: Optional[
        List[str]
    ] = "gs://broad-references-private/HybSelOligos/xgen_plus_spikein/white_album_exome_calling_regions.v1.interval_list",
    apply_recalibration_machine_mem_gb: Optional[int] = 60,
    is_small_callset: Optional[bool] = False,
):
    b = hb.Batch("VariantCallingOFTHEFUTURE")
    dbsnp_resource_vcf = (
        dbsnp_resource_vcf if dbsnp_resource_vcf is not None else dbsnp_vcf
    )
    dbsnp_resource_vcf_index = (
        dbsnp_resource_vcf_index
        if dbsnp_resource_vcf_index is not None
        else dbsnp_vcf_index
    )
    unboundedScatterCount = (
        unboundedScatterCount if unboundedScatterCount is not None else vcf_count
    )
    scatterCount = (
        scatterCount
        if scatterCount is not None
        else (
            unboundedScatterCount
            if (unboundedScatterCount and (unboundedScatterCount > 10))
            else 10
        )
    )
    combined_gvcf = b.read_input(combined_gvcf)
    ref_fasta = b.read_input_group(
        base=ref_fasta,
        dict=ref_fasta.replace(".fasta", "")
        .replace(".fna", "")
        .replace(".fa", "")
        + ".dict",
        fai=ref_fasta + ".fai",
    )
    ref_fasta_index = b.read_input(ref_fasta_index)
    ref_dict = b.read_input(ref_dict)
    dbsnp_vcf = b.read_input(dbsnp_vcf)
    dbsnp_vcf_index = b.read_input(dbsnp_vcf_index)
    eval_interval_list = b.read_input(eval_interval_list)
    hapmap_resource_vcf = b.read_input(hapmap_resource_vcf)
    hapmap_resource_vcf_index = b.read_input(hapmap_resource_vcf_index)
    omni_resource_vcf = b.read_input(omni_resource_vcf)
    omni_resource_vcf_index = b.read_input(omni_resource_vcf_index)
    one_thousand_genomes_resource_vcf = b.read_input(
        one_thousand_genomes_resource_vcf
    )
    one_thousand_genomes_resource_vcf_index = b.read_input(
        one_thousand_genomes_resource_vcf_index
    )
    mills_resource_vcf = b.read_input(mills_resource_vcf)
    mills_resource_vcf_index = b.read_input(mills_resource_vcf_index)
    axiomPoly_resource_vcf = b.read_input(axiomPoly_resource_vcf)
    axiomPoly_resource_vcf_index = b.read_input(axiomPoly_resource_vcf_index)
    dbsnp_resource_vcf = b.read_input(dbsnp_resource_vcf)
    dbsnp_resource_vcf_index = b.read_input(dbsnp_resource_vcf_index)
    unpadded_intervals = [
        b.read_input(inner_unpadded_intervals)
        for inner_unpadded_intervals in unpadded_intervals
    ]

    TabixBGzippedFile = add_TabixBGzippedFile_step(b, inp=combined_gvcf)

    GnarlyGenotyperOnVcf = []
    for idx in range(len(unpadded_intervals)):
        GnarlyGenotyperOnVcf.append(
            add_GnarlyGenotyperOnVcf_step(
                b,
                combined_gvcf=TabixBGzippedFile.out,
                interval=unpadded_intervals[idx],
                output_vcf_filename=(
                    ((callset_name + ".") + str(idx)) + ".vcf.gz"
                ),
                ref_fasta=ref_fasta,
                ref_fasta_index=ref_fasta_index,
                ref_dict=ref_dict,
                dbsnp_vcf=dbsnp_vcf,
                is_small_callset=is_small_callset,
            )
        )

    HardFilterAndMakeSitesOnlyVcf = []
    for idx in range(len(unpadded_intervals)):
        HardFilterAndMakeSitesOnlyVcf.append(
            add_HardFilterAndMakeSitesOnlyVcf_step(
                b,
                vcf=GnarlyGenotyperOnVcf[idx].output_vcf,
                excess_het_threshold=excess_het_threshold,
                variant_filtered_vcf_filename=(
                    ((callset_name + ".") + str(idx)) + ".variant_filtered.vcf.gz"
                ),
                sites_only_vcf_filename=(
                    ((callset_name + ".") + str(idx))
                    + ".sites_only.variant_filtered.vcf.gz"
                ),
                disk_size=medium_disk,
            )
        )
    SitesOnlyGatherVcf = add_SitesOnlyGatherVcf_step(
        b,
        input_vcfs=[j.sites_only_vcf for j in HardFilterAndMakeSitesOnlyVcf],
        output_vcf_name=(callset_name + ".sites_only.vcf.gz"),
        disk_size=medium_disk,
    )
    IndelsVariantRecalibrator = add_IndelsVariantRecalibrator_step(
        b,
        sites_only_variant_filtered_vcf=SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index=SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename=(callset_name + ".indels.recal"),
        tranches_filename=(callset_name + ".indels.tranches"),
        recalibration_tranche_values=indel_recalibration_tranche_values,
        recalibration_annotation_values=indel_recalibration_annotation_values,
        mills_resource_vcf=mills_resource_vcf,
        mills_resource_vcf_index=mills_resource_vcf_index,
        axiomPoly_resource_vcf=axiomPoly_resource_vcf,
        axiomPoly_resource_vcf_index=axiomPoly_resource_vcf_index,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        dbsnp_resource_vcf_index=dbsnp_resource_vcf_index,
        use_allele_specific_annotations=use_allele_specific_annotations,
        disk_size=small_disk,
    )
    SNPsVariantRecalibratorCreateModel = add_SNPsVariantRecalibratorCreateModel_step(
        b,
        sites_only_variant_filtered_vcf=SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index=SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename=(callset_name + ".snps.recal"),
        tranches_filename=(callset_name + ".snps.tranches"),
        recalibration_tranche_values=snp_recalibration_tranche_values,
        recalibration_annotation_values=snp_recalibration_annotation_values,
        downsampleFactor=SNP_VQSR_downsampleFactor,
        model_report_filename=(callset_name + ".snps.model.report"),
        hapmap_resource_vcf=hapmap_resource_vcf,
        hapmap_resource_vcf_index=hapmap_resource_vcf_index,
        omni_resource_vcf=omni_resource_vcf,
        omni_resource_vcf_index=omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index=one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        dbsnp_resource_vcf_index=dbsnp_resource_vcf_index,
        disk_size=small_disk,
        use_allele_specific_annotations=use_allele_specific_annotations,
    )

    SNPsVariantRecalibratorScattered = []
    for idx in range(len(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf)):
        SNPsVariantRecalibratorScattered.append(
            add_SNPsVariantRecalibratorScattered_step(
                b,
                sites_only_variant_filtered_vcf=HardFilterAndMakeSitesOnlyVcf.sites_only_vcf[
                    str(idx)
                ],
                sites_only_variant_filtered_vcf_index=HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index[
                    str(idx)
                ],
                recalibration_filename=(
                    ((callset_name + ".snps.") + str(idx)) + ".recal"
                ),
                tranches_filename=(
                    ((callset_name + ".snps.") + str(idx)) + ".tranches"
                ),
                recalibration_tranche_values=snp_recalibration_tranche_values,
                recalibration_annotation_values=snp_recalibration_annotation_values,
                model_report=SNPsVariantRecalibratorCreateModel.model_report,
                hapmap_resource_vcf=hapmap_resource_vcf,
                hapmap_resource_vcf_index=hapmap_resource_vcf_index,
                omni_resource_vcf=omni_resource_vcf,
                omni_resource_vcf_index=omni_resource_vcf_index,
                one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
                one_thousand_genomes_resource_vcf_index=one_thousand_genomes_resource_vcf_index,
                dbsnp_resource_vcf=dbsnp_resource_vcf,
                dbsnp_resource_vcf_index=dbsnp_resource_vcf_index,
                disk_size=small_disk,
                machine_mem_gb=apply_recalibration_machine_mem_gb,
                use_allele_specific_annotations=use_allele_specific_annotations,
            )
        )
    SNPGatherTranches = add_SNPGatherTranches_step(
        b,
        tranches=SNPsVariantRecalibratorScattered.tranches,
        output_filename=(callset_name + ".snps.gathered.tranches"),
        disk_size=small_disk,
    )

    ApplyRecalibration = []
    for idx in range(len(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf)):
        ApplyRecalibration.append(
            add_ApplyRecalibration_step(
                b,
                recalibrated_vcf_filename=(
                    ((callset_name + ".filtered.") + idx) + ".vcf.gz"
                ),
                input_vcf=HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
                input_vcf_index=HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[
                    idx
                ],
                indels_recalibration=IndelsVariantRecalibrator.recalibration,
                indels_recalibration_index=IndelsVariantRecalibrator.recalibration_index,
                indels_tranches=IndelsVariantRecalibrator.tranches,
                snps_recalibration=SNPsVariantRecalibratorScattered.recalibration[
                    idx
                ],
                snps_recalibration_index=SNPsVariantRecalibratorScattered.recalibration_index[
                    idx
                ],
                snps_tranches=SNPGatherTranches.out_tranches,
                indel_filter_level=indel_filter_level,
                snp_filter_level=snp_filter_level,
                disk_size=medium_disk,
                use_allele_specific_annotations=use_allele_specific_annotations,
            )
        )
    CollectMetricsSharded = None
    if not is_small_callset:
        CollectMetricsSharded = add_CollectMetricsSharded_step(
            b,
            input_vcf=ApplyRecalibration.recalibrated_vcf,
            input_vcf_index=ApplyRecalibration.recalibrated_vcf_index,
            metrics_filename_prefix=((callset_name + ".") + idx),
            dbsnp_vcf=dbsnp_vcf,
            dbsnp_vcf_index=dbsnp_vcf_index,
            interval_list=eval_interval_list,
            ref_dict=ref_dict,
            disk_size=small_disk,
        )

    FinalGatherVcf = None
    if is_small_callset:
        FinalGatherVcf = add_FinalGatherVcf_step(
            b,
            input_vcfs=ApplyRecalibration.recalibrated_vcf,
            output_vcf_name=(callset_name + ".vcf.gz"),
            disk_size=huge_disk,
        )

    CollectMetricsOnFullVcf = None
    if is_small_callset:
        CollectMetricsOnFullVcf = add_CollectMetricsOnFullVcf_step(
            b,
            input_vcf=FinalGatherVcf.output_vcf,
            input_vcf_index=FinalGatherVcf.output_vcf_index,
            metrics_filename_prefix=callset_name,
            dbsnp_vcf=dbsnp_vcf,
            dbsnp_vcf_index=dbsnp_vcf_index,
            interval_list=eval_interval_list,
            ref_dict=ref_dict,
            disk_size=huge_disk,
        )

    return b


def add_TabixBGzippedFile_step(
    b, inp, preset="vcf", container="quay.io/biocontainers/htslib:1.9--ha228f0b_7"
):
    j = b.new_job("TabixBGzippedFile")
    j.command(f"mv {inp} .")
    j.declare_resource_group(out={"tbi": "{root}.tbi", "base": "{root}"})
    j.image(container)

    command_args = []
    if preset is not None:
        command_args.append(f"--preset {preset}")
    if inp is not None:
        command_args.append(inp)
    nl = " \\\n  "
    command = f"""
tabix {"".join(nl + a for a in command_args)} 
    """

    j.command('ln "{value}" {dest}'.format(value=inp, dest=j.out))
    j.command(
        'ln "{value}" {dest}'.format(value=f"{inp}.tbi", dest=f"{j.out}.tbi")
    )

    return j


def add_GnarlyGenotyperOnVcf_step(
    b,
    combined_gvcf,
    interval,
    output_vcf_filename,
    ref_fasta,
    ref_fasta_index,
    ref_dict,
    dbsnp_vcf,
    gatk_docker="gcr.io/broad-dsde-methods/gnarly_genotyper:hail_ukbb_300K",
    is_small_callset=False,
    disk_size=None,
    container="gcr.io/broad-dsde-methods/gnarly_genotyper:hail_ukbb_300K",
):
    j = b.new_job("GnarlyGenotyperOnVcf")
    disk_size = (
        disk_size if disk_size is not None else (40 if is_small_callset else 80)
    )
    j.declare_resource_group(output_vcf={"tbi": "{root}.tbi", "base": "{root}"})
    j.image(container)
    j.memory(f"24.214398G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

    j.command(
        f"""set -e

    gatk --java-options -Xms8g \\
      GnarlyGenotyper \\
      -R {ref_fasta} \\
      -O {output_vcf_filename} \\
      -D {dbsnp_vcf} \\
      --only-output-calls-starting-in-intervals \\
      --keep-all-sites \\
      -V {combined_gvcf} \\
      -L {interval}"""
    )

    j.command(
        'ln "{value}" {dest}'.format(value=output_vcf_filename, dest=j.output_vcf)
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{output_vcf_filename}.tbi", dest=f"{j.output_vcf}.tbi"
        )
    )

    return j


def add_HardFilterAndMakeSitesOnlyVcf_step(
    b,
    vcf,
    excess_het_threshold,
    variant_filtered_vcf_filename,
    sites_only_vcf_filename,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("HardFilterAndMakeSitesOnlyVcf")
    j.image(container)
    j.memory(f"3.49246125G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms3g \\
      VariantFiltration \\
      --filter-expression 'ExcessHet > {excess_het_threshold}' \\
      --filter-name ExcessHet \\
      -O {variant_filtered_vcf_filename} \\
      -V {vcf.base}

    gatk --java-options -Xms3g \\
      MakeSitesOnlyVcf \\
      -I {variant_filtered_vcf_filename} \\
      -O {sites_only_vcf_filename}"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{variant_filtered_vcf_filename}", dest=j.variant_filtered_vcf
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{variant_filtered_vcf_filename}.tbi",
            dest=j.variant_filtered_vcf_index,
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{sites_only_vcf_filename}", dest=j.sites_only_vcf
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{sites_only_vcf_filename}.tbi", dest=j.sites_only_vcf_index
        )
    )

    return j


def add_SitesOnlyGatherVcf_step(
    b,
    input_vcfs,
    output_vcf_name,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("SitesOnlyGatherVcf")
    j.image(container)
    j.memory(f"6.519261G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

    j.command(
        f"""set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options -Xms6g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      --input {input_vcfs} \\
      --output {output_vcf_name}

    tabix {output_vcf_name}"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{output_vcf_name}", dest=j.output_vcf
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{output_vcf_name}.tbi", dest=j.output_vcf_index
        )
    )

    return j


def add_IndelsVariantRecalibrator_step(
    b,
    recalibration_filename,
    tranches_filename,
    recalibration_tranche_values,
    recalibration_annotation_values,
    sites_only_variant_filtered_vcf,
    sites_only_variant_filtered_vcf_index,
    mills_resource_vcf,
    axiomPoly_resource_vcf,
    dbsnp_resource_vcf,
    mills_resource_vcf_index,
    axiomPoly_resource_vcf_index,
    dbsnp_resource_vcf_index,
    use_allele_specific_annotations,
    disk_size,
    model_report=None,
    max_gaussians=4,
    model_report_arg=None,
    container="ubuntu:latest",
):
    j = b.new_job("IndelsVariantRecalibrator")
    model_report_arg = (
        model_report_arg
        if model_report_arg is not None
        else (
            "--input-model $MODEL_REPORT --output-tranches-for-scatter"
            if model_report is not None
            else ""
        )
    )
    j.image(container)
    j.memory(f"104.0G")
    j.storage(f"disk_sizeG")

    j.command(
        f"""set -euo pipefail

    MODEL_REPORT={model_report.base}

    gatk --java-options -Xms100g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf.base} \\
      -O {recalibration_filename} \\
      --tranches-file {tranches_filename} \\
      --trust-all-polymorphic \\
      -tranche {recalibration_tranche_values} \\
      -an {recalibration_annotation_values} \\
      -mode INDEL \\
      {use_allele_specific_annotations} \\
      {model_report_arg} \\
      --max-gaussians {max_gaussians} \\
      -resource:mills,known=false,training=true,truth=true,prior=12 {mills_resource_vcf.base} \\
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {axiomPoly_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 {dbsnp_resource_vcf.base}"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{recalibration_filename}", dest=j.recalibration
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{recalibration_filename}.idx", dest=j.recalibration_index
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{tranches_filename}", dest=j.tranches
        )
    )

    return j


def add_SNPsVariantRecalibratorCreateModel_step(
    b,
    recalibration_filename,
    tranches_filename,
    downsampleFactor,
    model_report_filename,
    recalibration_tranche_values,
    recalibration_annotation_values,
    sites_only_variant_filtered_vcf,
    sites_only_variant_filtered_vcf_index,
    hapmap_resource_vcf,
    omni_resource_vcf,
    one_thousand_genomes_resource_vcf,
    dbsnp_resource_vcf,
    hapmap_resource_vcf_index,
    omni_resource_vcf_index,
    one_thousand_genomes_resource_vcf_index,
    dbsnp_resource_vcf_index,
    disk_size,
    use_allele_specific_annotations,
    max_gaussians=6,
    java_mem=100,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.4.1",
    container="us.gcr.io/broad-gatk/gatk:4.1.4.1",
):
    j = b.new_job("SNPsVariantRecalibratorCreateModel")
    j.image(container)
    j.memory(f"104.0G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{java_mem}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf.base} \\
      -O {recalibration_filename} \\
      --tranches-file {tranches_filename} \\
      --trust-all-polymorphic \\
      -tranche {recalibration_tranche_values} \\
      -an {recalibration_annotation_values} \\
      -mode SNP \\
      {use_allele_specific_annotations} \\
      --sample-every-Nth-variant {downsampleFactor} \\
      --output-model {model_report_filename} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf}"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{model_report_filename}", dest=j.model_report
        )
    )

    return j


def add_SNPsVariantRecalibratorScattered_step(
    b,
    recalibration_filename,
    tranches_filename,
    recalibration_tranche_values,
    recalibration_annotation_values,
    sites_only_variant_filtered_vcf,
    sites_only_variant_filtered_vcf_index,
    hapmap_resource_vcf,
    omni_resource_vcf,
    one_thousand_genomes_resource_vcf,
    dbsnp_resource_vcf,
    hapmap_resource_vcf_index,
    omni_resource_vcf_index,
    one_thousand_genomes_resource_vcf_index,
    dbsnp_resource_vcf_index,
    disk_size,
    use_allele_specific_annotations,
    model_report=None,
    max_gaussians=6,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    machine_mem_gb=None,
    auto_mem=None,
    machine_mem=None,
    model_report_arg=None,
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("SNPsVariantRecalibratorScattered")
    auto_mem = (
        auto_mem
        if auto_mem is not None
        else math.ceil(
            (
                2
                * (
                    None
                    * (
                        (
                            (
                                (
                                    (
                                        0
                                        + os.stat(
                                            sites_only_variant_filtered_vcf
                                        ).st_size
                                        / 1000
                                    )
                                    + os.stat(hapmap_resource_vcf).st_size / 1000
                                )
                                + os.stat(omni_resource_vcf).st_size / 1000
                            )
                            + os.stat(one_thousand_genomes_resource_vcf).st_size
                            / 1000
                        )
                        + os.stat(dbsnp_resource_vcf).st_size / 1000
                    )
                )
            )
        )
    )
    machine_mem = (
        machine_mem
        if machine_mem is not None
        else [
            a
            for a in [machine_mem_gb, (7 if (auto_mem < 7) else auto_mem)]
            if a is not None
        ]
    )
    model_report_arg = (
        model_report_arg
        if model_report_arg is not None
        else (
            "--input-model $MODEL_REPORT --output-tranches-for-scatter"
            if model_report is not None
            else ""
        )
    )
    j.image(container)
    j.memory(f"{machine_mem} GiBG")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

    j.command(
        f"""set -euo pipefail

    MODEL_REPORT={model_report.base}

    gatk --java-options -Xms{{(machine_mem - 1)}}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf.base} \\
      -O {recalibration_filename} \\
      --tranches-file {tranches_filename} \\
      --trust-all-polymorphic \\
      -tranche {recalibration_tranche_values} \\
      -an {recalibration_annotation_values} \\
      -mode SNP \\
      {use_allele_specific_annotations} \\
       {model_report_arg} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base}"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{recalibration_filename}", dest=j.recalibration
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{recalibration_filename}.idx", dest=j.recalibration_index
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{tranches_filename}", dest=j.tranches
        )
    )

    return j


def add_SNPGatherTranches_step(
    b,
    tranches,
    output_filename,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("SNPGatherTranches")
    j.image(container)
    j.memory(f"6.519261G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

    j.command(
        f"""set -euo pipefail

    tranches_fofn={JANIS: write_lines([inputs.tranches])}

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
      --output {output_filename}"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{output_filename}", dest=j.out_tranches
        )
    )

    return j


def add_ApplyRecalibration_step(
    b,
    recalibrated_vcf_filename,
    input_vcf,
    input_vcf_index,
    indels_recalibration,
    indels_recalibration_index,
    indels_tranches,
    snps_recalibration,
    snps_recalibration_index,
    snps_tranches,
    indel_filter_level,
    snp_filter_level,
    use_allele_specific_annotations,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("ApplyRecalibration")
    j.image(container)
    j.memory(f"7.0G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O tmp.indel.recalibrated.vcf \\
      -V {input_vcf.base} \\
      --recal-file {indels_recalibration.base} \\
      --tranches-file {indels_tranches.base} \\
      --truth-sensitivity-filter-level {indel_filter_level} \\
      --create-output-variant-index true \\
      -mode INDEL {use_allele_specific_annotations} \\


    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O {recalibrated_vcf_filename} \\
      -V tmp.indel.recalibrated.vcf \\
      --recal-file {snps_recalibration.base} \\
      --tranches-file {snps_tranches.base} \\
      --truth-sensitivity-filter-level {snp_filter_level} \\
      --create-output-variant-index true \\
      -mode SNP {use_allele_specific_annotations} \\
"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{recalibrated_vcf_filename}", dest=j.recalibrated_vcf
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{recalibrated_vcf_filename}.tbi",
            dest=j.recalibrated_vcf_index,
        )
    )

    return j


def add_CollectMetricsSharded_step(
    b,
    input_vcf,
    input_vcf_index,
    metrics_filename_prefix,
    dbsnp_vcf,
    dbsnp_vcf_index,
    interval_list,
    ref_dict,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("CollectMetricsSharded")
    j.image(container)
    j.memory(f"6.9849225G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      CollectVariantCallingMetrics \\
      --INPUT {input_vcf.base} \\
      --DBSNP {dbsnp_vcf.base} \\
      --SEQUENCE_DICTIONARY {ref_dict.base} \\
      --OUTPUT {metrics_filename_prefix} \\
      --THREAD_COUNT 8 \\
      --TARGET_INTERVALS {interval_list.base}"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{metrics_filename_prefix}.variant_calling_detail_metrics",
            dest=j.detail_metrics_file,
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{metrics_filename_prefix}.variant_calling_summary_metrics",
            dest=j.summary_metrics_file,
        )
    )

    return j


def add_FinalGatherVcf_step(
    b,
    input_vcfs,
    output_vcf_name,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("FinalGatherVcf")
    j.image(container)
    j.memory(f"6.519261G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

    j.command(
        f"""set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options -Xms6g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      --input {input_vcfs} \\
      --output {output_vcf_name}

    tabix {output_vcf_name}"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{output_vcf_name}", dest=j.output_vcf
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{output_vcf_name}.tbi", dest=j.output_vcf_index
        )
    )

    return j


def add_CollectMetricsOnFullVcf_step(
    b,
    input_vcf,
    input_vcf_index,
    metrics_filename_prefix,
    dbsnp_vcf,
    dbsnp_vcf_index,
    interval_list,
    ref_dict,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("CollectMetricsOnFullVcf")
    j.image(container)
    j.memory(f"6.9849225G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      CollectVariantCallingMetrics \\
      --INPUT {input_vcf.base} \\
      --DBSNP {dbsnp_vcf.base} \\
      --SEQUENCE_DICTIONARY {ref_dict.base} \\
      --OUTPUT {metrics_filename_prefix} \\
      --THREAD_COUNT 8 \\
      --TARGET_INTERVALS {interval_list.base}"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{metrics_filename_prefix}.variant_calling_detail_metrics",
            dest=j.detail_metrics_file,
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{metrics_filename_prefix}.variant_calling_summary_metrics",
            dest=j.summary_metrics_file,
        )
    )

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


@click.command()
@click.option("--combined_gvcf", "combined_gvcf", type=str, required=True)
@click.option("--callset_name", "callset_name", type=str, required=True)
@click.option(
    "--unpadded_intervals_file",
    "unpadded_intervals_file",
    type=str,
    default="gs://broad-references-private/HybSelOligos/xgen_plus_spikein/white_album_exome_calling_regions.v1.interval_list",
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
    default="gs://gcp-public-data--broad-references/hg38/v0/exome_evaluation_regions.v1.interval_list",
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
    "--axiomPoly_resource_vcf",
    "axiomPoly_resource_vcf",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
)
@click.option(
    "--axiomPoly_resource_vcf_index",
    "axiomPoly_resource_vcf_index",
    type=str,
    default="gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi",
)
@click.option("--dbsnp_resource_vcf", "dbsnp_resource_vcf", type=str)
@click.option("--dbsnp_resource_vcf_index", "dbsnp_resource_vcf_index", type=str)
@click.option(
    "--excess_het_threshold", "excess_het_threshold", type=float, default=54.69
)
@click.option("--snp_filter_level", "snp_filter_level", type=float, default=99.7)
@click.option(
    "--indel_filter_level", "indel_filter_level", type=float, default=99.0
)
@click.option(
    "--SNP_VQSR_downsampleFactor",
    "SNP_VQSR_downsampleFactor",
    type=int,
    default=75,
)
@click.option(
    "--indel_VQSR_downsampleFactor",
    "indel_VQSR_downsampleFactor",
    type=int,
    default=10,
)
@click.option(
    "--use_allele_specific_annotations",
    "use_allele_specific_annotations",
    is_flag=True,
    default=True,
)
@click.option("--vcf_count", "vcf_count", type=int, default=30)
@click.option("--unboundedScatterCount", "unboundedScatterCount", type=int)
@click.option("--scatterCount", "scatterCount", type=int)
@click.option(
    "--unpadded_intervals",
    "unpadded_intervals",
    multiple=True,
    type=str,
    default="gs://broad-references-private/HybSelOligos/xgen_plus_spikein/white_album_exome_calling_regions.v1.interval_list",
)
@click.option(
    "--apply_recalibration_machine_mem_gb",
    "apply_recalibration_machine_mem_gb",
    type=int,
    default=60,
)
@click.option(
    "--is_small_callset", "is_small_callset", is_flag=True, default=False
)
def main_from_click(*args, **kwargs):
    return main(*args, **kwargs)


if __name__ == "__main__":
    main_from_click()
