import math, os, tempfile
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
    backend = hb.ServiceBackend(billing_project="test")
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
    axiomPoly_resource_vcf = b.read_input_group(
        base=axiomPoly_resource_vcf, 
        tbi=axiomPoly_resource_vcf_index)
    dbsnp_resource_vcf = b.read_input_group(
        base=dbsnp_resource_vcf, 
        tbi=dbsnp_resource_vcf_index)

    unpadded_intervals = [
        b.read_input(inner_unpadded_intervals)
        for inner_unpadded_intervals in unpadded_intervals
    ]
    unpadded_intervals = unpadded_intervals[0]

    GnarlyGenotyperOnVcf = []
    for idx in range(len(unpadded_intervals)):
        GnarlyGenotyperOnVcf.append(
            add_GnarlyGenotyperOnVcf_step(
                b,
                combined_gvcf=combined_gvcf,
                interval=unpadded_intervals[idx],
                output_vcf_filename=(
                    ((callset_name + ".") + str(idx)) + ".vcf.gz"
                ),
                ref_fasta=ref_fasta,
                dbsnp_vcf=dbsnp_vcf,
                is_small_callset=is_small_callset,
            )
        )

    HardFilterAndMakeSitesOnlyVcf = []
    for idx in range(len(unpadded_intervals)):
        HardFilterAndMakeSitesOnlyVcf.append(
            add_HardFilterAndMakeSitesOnlyVcf_step(
                b,
                input_vcfs=[j.output_vcf for j in GnarlyGenotyperOnVcf],
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
        recalibration_filename=(callset_name + ".indels.recal"),
        tranches_filename=(callset_name + ".indels.tranches"),
        recalibration_tranche_values=indel_recalibration_tranche_values,
        recalibration_annotation_values=indel_recalibration_annotation_values,
        mills_resource_vcf=mills_resource_vcf,
        axiomPoly_resource_vcf=axiomPoly_resource_vcf,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        use_allele_specific_annotations=use_allele_specific_annotations,
        disk_size=small_disk,
    )
    SNPsVariantRecalibratorCreateModel = add_SNPsVariantRecalibratorCreateModel_step(
        b,
        sites_only_variant_filtered_vcf=SitesOnlyGatherVcf.output_vcf,
        recalibration_filename=(callset_name + ".snps.recal"),
        tranches_filename=(callset_name + ".snps.tranches"),
        recalibration_tranche_values=snp_recalibration_tranche_values,
        recalibration_annotation_values=snp_recalibration_annotation_values,
        downsampleFactor=SNP_VQSR_downsampleFactor,
        model_report_filename=(callset_name + ".snps.model.report"),
        hapmap_resource_vcf=hapmap_resource_vcf,
        omni_resource_vcf=omni_resource_vcf,
        one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        disk_size=small_disk,
        use_allele_specific_annotations=use_allele_specific_annotations,
    )

    SNPsVariantRecalibratorScattered = []
    for j in HardFilterAndMakeSitesOnlyVcf:
        SNPsVariantRecalibratorScattered.append(
            add_SNPsVariantRecalibratorScattered_step(
                b,
                sites_only_variant_filtered_vcf=j.sites_only_vcf,
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
                omni_resource_vcf=omni_resource_vcf,
                one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
                dbsnp_resource_vcf=dbsnp_resource_vcf,
                disk_size=small_disk,
                machine_mem_gb=apply_recalibration_machine_mem_gb,
                is_small_callset=is_small_callset,
                use_allele_specific_annotations=use_allele_specific_annotations,
            )
        )
    SNPGatherTranches = add_SNPGatherTranches_step(
        b,
        tranches=[j.tranches for j in SNPsVariantRecalibratorScattered],
        output_filename=(callset_name + ".snps.gathered.tranches"),
        disk_size=small_disk,
    )

    ApplyRecalibration = []
    for idx, j in enumerate(HardFilterAndMakeSitesOnlyVcf):
        ApplyRecalibration.append(
            add_ApplyRecalibration_step(
                b,
                recalibrated_vcf_filename=(
                    ((callset_name + ".filtered.") + str(idx)) + ".vcf.gz"
                ),
                input_vcf=j.variant_filtered_vcf,
                indels_recalibration=IndelsVariantRecalibrator.recalibration,
                indels_tranches=IndelsVariantRecalibrator.tranches,
                snps_recalibration=SNPsVariantRecalibratorScattered[idx].recalibration,
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
            metrics_filename_prefix=((callset_name + ".") + idx),
            dbsnp_vcf=dbsnp_vcf,
            interval_list=eval_interval_list,
            ref_dict=ref_fasta.dict,
            disk_size=small_disk,
        )

    FinalGatherVcf = None
    if is_small_callset:
        FinalGatherVcf = add_FinalGatherVcf_step(
            b,
            input_vcfs=[j.recalibrated_vcf for j in ApplyRecalibration],
            output_vcf_name=(callset_name + ".vcf.gz"),
            disk_size=huge_disk,
        )

    CollectMetricsOnFullVcf = None
    if is_small_callset:
        CollectMetricsOnFullVcf = add_CollectMetricsOnFullVcf_step(
            b,
            input_vcf=FinalGatherVcf.output_vcf,
            metrics_filename_prefix=callset_name,
            dbsnp_vcf=dbsnp_vcf,
            interval_list=eval_interval_list,
            ref_dict=ref_fasta.dict,
            disk_size=huge_disk,
        )

    b.run(dry_run=True)


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
    container="gcr.io/broad-dsde-methods/gnarly_genotyper:hail_ukbb_300K",
):
    j = b.new_job("GnarlyGenotyperOnVcf")
    disk_size = (
        disk_size if disk_size is not None else (40 if is_small_callset else 80)
    )
    j.declare_resource_group(output_vcf={"index": "{root}.tbi", "base": "{root}"})
    j.image(container)
    j.memory(f"24.214398G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

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
      -L {interval}"""
    )
    b.write_output(j.output_vcf, output_vcf_filename)
    return j


def add_HardFilterAndMakeSitesOnlyVcf_step(
    b,
    input_vcfs,
    excess_het_threshold,
    variant_filtered_vcf_filename,
    sites_only_vcf_filename,
    disk_size,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("HardFilterAndMakeSitesOnlyVcf")
    j.declare_resource_group(variant_filtered_vcf={"base": "{root}", "index": "{root}.tbi"})
    j.declare_resource_group(sites_only_vcf={"base": "{root}", "index": "{root}.tbi"})
    j.image(container)
    j.memory(f"3.49246125G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')

    input_cmdl = " ".join([v.base for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms3g \\
      VariantFiltration \\
      --filter-expression 'ExcessHet > {excess_het_threshold}' \\
      --filter-name ExcessHet \\
      -O {j.variant_filtered_vcf.base} \\
      -V {input_cmdl}

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
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("SitesOnlyGatherVcf")
    j.image(container)
    j.memory(f"6.519261G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')
    j.declare_resource_group(output_vcf={"base": "{root}", "index": "{root}.tbi"})

    input_cmdl = " ".join([v.base for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options -Xms6g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      --input {input_cmdl} \\
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
    axiomPoly_resource_vcf,
    dbsnp_resource_vcf,
    use_allele_specific_annotations,
    disk_size,
    max_gaussians=4,
    container="ubuntu:latest",
):
    j = b.new_job("IndelsVariantRecalibrator")
    j.image(container)
    j.memory(f"104.0G")
    j.storage(f"disk_sizeG")
    j.declare_resource_group(recalibration={"base": "{root}", "index": "{root}.idx"})

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms100g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf} \\
      -O {j.recalibration.base} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      -tranche {recalibration_tranche_values} \\
      -an {recalibration_annotation_values} \\
      -mode INDEL \\
      {use_allele_specific_annotations} \\
      --max-gaussians {max_gaussians} \\
      -resource:mills,known=false,training=true,truth=true,prior=12 {mills_resource_vcf.base} \\
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {axiomPoly_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 {dbsnp_resource_vcf.base}"""
    )
    b.write_output(j.recalibration, recalibration_filename)
    b.write_output(j.tranches, tranches_filename)
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
    hapmap_resource_vcf,
    omni_resource_vcf,
    one_thousand_genomes_resource_vcf,
    dbsnp_resource_vcf,
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
    j.declare_resource_group(recalibration={"base": "{root}", "index": "{root}.idx"})

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{java_mem}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf.base} \\
      -O {j.recalibration.base} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      -tranche {recalibration_tranche_values} \\
      -an {recalibration_annotation_values} \\
      -mode SNP \\
      {use_allele_specific_annotations} \\
      --sample-every-Nth-variant {downsampleFactor} \\
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
    is_small_callset,
    model_report=None,
    max_gaussians=6,
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    machine_mem_gb=None,
    machine_mem=None,
    model_report_arg=None,
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("SNPsVariantRecalibratorScattered")
    machine_mem = (
        machine_mem
        if machine_mem is not None
        else [
            a
            for a in [machine_mem_gb, (30 if is_small_callset else 60)]
            if a is not None
        ]
    )
    model_report_arg = (
        model_report_arg
        if model_report_arg is not None
        else (
            (("--input-model " + model_report) + "--output-tranches-for-scatter")
            if model_report is not None
            else ""
        )
    )
    j.image(container)
    j.memory(f"{machine_mem} GiBG")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')
    j.declare_resource_group(recalibration={"base": "{root}", "index": "{root}.index"})

    j.command(
        f"""set -euo pipefail

    MODEL_REPORT={model_report}

    gatk --java-options -Xms{{(machine_mem - 1)}}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf.base} \\
      -O {j.recalibration.base} \\
      --tranches-file {j.tranches} \\
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
    b.write_output(j.recalibration, recalibration_filename)
    b.write_output(j.tranches, tranches_filename)
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

    tranches_file = tempfile.NamedTemporaryFile(suffix="_tranches.list")
    with open(tranches_file.name, 'w') as out:
        for tranch in tranches:
            out.write(f'{tranch}\n')

    j.command(
        f"""set -euo pipefail

    tranches_fofn={tranches_file.name}

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
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("ApplyRecalibration")
    j.image(container)
    j.memory(f"7.0G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')
    j.declare_resource_group(recalibrated_vcf={"base": "{root}", "index": "{root}.tbi"})

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O tmp.indel.recalibrated.vcf \\
      -V {input_vcf.base} \\
      --recal-file {indels_recalibration.base} \\
      --tranches-file {indels_tranches} \\
      --truth-sensitivity-filter-level {indel_filter_level} \\
      --create-output-variant-index true \\
      -mode INDEL {use_allele_specific_annotations} \\

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O {j.recalibrated_vcf.base} \\
      -V tmp.indel.recalibrated.vcf \\
      --recal-file {snps_recalibration.base} \\
      --tranches-file {snps_tranches} \\
      --truth-sensitivity-filter-level {snp_filter_level} \\
      --create-output-variant-index true \\
      -mode SNP {use_allele_specific_annotations} \\
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
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("CollectMetricsSharded")
    j.image(container)
    j.memory(f"6.9849225G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')
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
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("FinalGatherVcf")
    j.image(container)
    j.memory(f"6.519261G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')
    j.declare_resource_group(output_vcf = {"base": "{root}", "index": "{root}.tbi"})

    input_cmdl = " ".join([v.base for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options -Xms6g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      --input {input_cmdl} \\
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
    gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
):
    j = b.new_job("CollectMetricsOnFullVcf")
    j.image(container)
    j.memory(f"6.9849225G")
    j.storage(f'{(("local-disk " + str(disk_size)) + " HDD")}G')
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
