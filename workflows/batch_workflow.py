"""
Hail Batch workflow to perform joint calling, sample QC, and variant QC with VQSR and 
random forest methods on a WGS germline callset.
1. For the sample QC and the random forest variant QC, mostly re-implementation of 
https://github.com/broadinstitute/gnomad_qc
2. For VQSR, compilation from the following 2 WDL workflows:
a. hail-ukbb-200k-callset/GenotypeAndFilter.AS.wdl
b. The Broad VQSR workflow:
   https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/dna_seq/germline/joint_genotyping/JointGenotyping.wdl
   documented here:
   https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
Translated from WDL with a help of Janis:
https://github.com/PMCC-BioinformaticsCore/janis

Output
* The output of sample QC is a CSV file:
`<output_bucket>/meta.tsv`
* The output of VQSR is a VCF file <output_bucket>/<callset>-recalibrated.vcf.gz,
as well as a QC file <output_bucket>/<callset>-eval.txt
and R scripts to plot VQSR models: <output_bucket>/plot-snps-recal.Rscript
and <output_bucket>/plot-indels-recal.Rscript

The workflow is parametrised by the dataset name, batch names and the output version.
$ python scripts/batch_workflow.py --callset fewgenomes \
     --version v0 --batch 0 --batch 1 --billing-project test --skip-qc
Will read the inputs from:
gs://cpg-fewgenomes-test/gvcf/batch0/*.g.vcf.gz
And write into the following output bucket:
gs://cpg-fewgenomes-temporary/joint_vcf/v0/work/combiner
"""

import os
import uuid
from os.path import join
from typing import List
import logging
import click
import pandas as pd
import subprocess
import hailtop.batch as hb
from hailtop.batch.job import Job
from analysis_runner import dataproc

from joint_calling import utils
from joint_calling.vqsr import make_vqsr_jobs

logger = logging.getLogger('joint-calling')
logger.setLevel('INFO')


@click.command()
@click.option('--callset', 'callset_name', type=str, required=True)
@click.option('--version', 'callset_version', type=str, required=True)
@click.option('--batch', 'callset_batches', type=str, multiple=True, required=True)
@click.option(
    '--from',
    'input_bucket_suffix',
    type=click.Choice(['main', 'test']),
    default='test',
    help='The bucket type to read from (default: test)',
)
@click.option(
    '--to',
    'output_bucket_suffix',
    type=click.Choice(['analysis', 'temporary']),
    default='temporary',
    help='The bucket type to write to (default: temporary)',
)
@click.option('--skip-qc', 'skip_qc', is_flag=True)
@click.option('--keep-scratch', 'keep_scratch', is_flag=True)
@click.option(
    '--reuse-scratch-run-id',
    'reuse_scratch_run_id',
    help='Run ID to reuse scratch from',
)
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option(
    '--billing-project',
    'billing_project',
    type=str,
    default=os.getenv('HAIL_BILLING_PROJECT'),
)
@click.option(
    '--snp_recalibration_tranche_values',
    'snp_recalibration_tranche_values',
    multiple=True,
    type=float,
    default=[
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
    ],
)
@click.option(
    '--snp_recalibration_annotation_values',
    'snp_recalibration_annotation_values',
    multiple=True,
    type=str,
    default=[
        'AS_QD',
        'AS_MQRankSum',
        'AS_ReadPosRankSum',
        'AS_FS',
        'AS_SOR',
        'AS_MQ',
    ],
)
@click.option(
    '--indel_recalibration_tranche_values',
    'indel_recalibration_tranche_values',
    multiple=True,
    type=float,
    default=[
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
    ],
)
@click.option(
    '--indel_recalibration_annotation_values',
    'indel_recalibration_annotation_values',
    multiple=True,
    type=str,
    default=['AS_FS', 'AS_SOR', 'AS_ReadPosRankSum', 'AS_MQRankSum', 'AS_QD'],
)
@click.option(
    '--excess_het_threshold', 'excess_het_threshold', type=float, default=54.69
)
@click.option('--snp_filter_level', 'snp_filter_level', type=float, default=99.7)
@click.option('--indel_filter_level', 'indel_filter_level', type=float, default=99.0)
@click.option(
    '--skip_allele_specific_annotations',
    'skip_allele_specific_annotations',
    is_flag=True,
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,too-many-statements
    callset_name: str,
    callset_version: str,
    callset_batches: List[str],
    input_bucket_suffix: str,
    output_bucket_suffix: str,
    skip_qc: bool,
    keep_scratch: bool,
    reuse_scratch_run_id: str,
    dry_run: bool,
    billing_project: str,
    snp_recalibration_tranche_values: List[float],
    snp_recalibration_annotation_values: List[str],
    indel_recalibration_tranche_values: List[float],
    indel_recalibration_annotation_values: List[str],
    excess_het_threshold: float,
    snp_filter_level: float,
    indel_filter_level: float,
    skip_allele_specific_annotations: bool,
):
    """
    Drive a Hail Batch workflow that creates and submits jobs. A job usually runs
    either a Hail Query script from the scripts folder in this repo using a Dataproc
    cluster; or a GATK command using the GATK or Gnarly image.
    """
    if not dry_run:
        if not billing_project:
            raise click.BadParameter(
                '--billing-project has to be specified (unless --dry-run is set)'
            )

    input_buckets = []
    for cb in callset_batches:
        cb = f'batch{cb}' if not cb.startswith('batch') else cb
        input_buckets.append(f'gs://cpg-{callset_name}-{input_bucket_suffix}/gvcf/{cb}')

    base_bucket = f'gs://cpg-{callset_name}-{output_bucket_suffix}/joint_vcf'
    output_bucket = join(base_bucket, f'{callset_version}/work')
    hail_bucket = join(output_bucket, 'hail')

    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=hail_bucket.replace('gs://', ''),
    )
    b = hb.Batch('Joint Calling', backend=backend)

    # TODO: merge with existing data

    samples_path = join(output_bucket, 'samples.csv')
    if not utils.file_exists(samples_path):
        samples_df = utils.find_inputs(input_buckets, skip_qc=skip_qc)
    else:
        samples_df = pd.read_csv(samples_path, sep='\t').set_index('s', drop=False)
    samples_df = samples_df[pd.notnull(samples_df.s)]

    gvcfs = [
        b.read_input_group(**{'g.vcf.gz': gvcf, 'g.vcf.gz.tbi': gvcf + '.tbi'})
        for gvcf in list(samples_df.gvcf)
    ]

    is_small_callset = len(gvcfs) < 1000
    # 1. For small callsets, we don't apply the ExcessHet filtering.
    # 2. For small callsets, we gather the VCF shards and collect QC metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics
    # collected from them.
    is_huge_callset = len(gvcfs) >= 100000
    # For huge callsets, we allocate more memory for the SNPs Create Model step

    # pylint: disable=unused-variable
    noalt_regions = b.read_input('gs://cpg-reference/hg38/v0/noalt.bed')

    # reblocked_gvcf_paths = [
    #     join(
    #         hail_bucket,
    #         'batch',
    #         reuse_scratch_run_id,
    #         str(job_num),
    #         'output_gvcf.g.vcf.gz',
    #     )
    #     if reuse_scratch_run_id
    #     else None
    #     for job_num in (1, 1 + len(gvcfs))
    # ]
    # reblocked_gvcfs = [
    #     b.read_input_group(
    #         **{
    #             'vcf.gz': output_gvcf_path,
    #             'vcf.gz.tbi': output_gvcf_path + '.tbi',
    #         }
    #     )
    #     if output_gvcf_path and utils.file_exists(output_gvcf_path)
    #     else add_reblock_gvcfs_step(b, input_gvcf).output_gvcf
    #     for output_gvcf_path, input_gvcf in zip(reblocked_gvcf_paths, gvcfs)
    # ]

    combiner_bucket = os.path.join(output_bucket, 'combiner')
    combiner_gvcf_bucket = os.path.join(output_bucket, 'combiner', 'gvcfs')
    # subset_gvcf_jobs = [
    #     add_subset_noalt_step(
    #         b,
    #         input_gvcf=gvcf,
    #         output_gvcf_path=join(combiner_gvcf_bucket, f'{sample}.g.vcf.gz'),
    #         noalt_regions=noalt_regions,
    #     )
    #     for sample, gvcf in zip(list(samples_df.s), reblocked_gvcfs)
    # ]
    for sn in samples_df.s:
        samples_df.loc[sn, ['gvcf']] = join(combiner_gvcf_bucket, sn + '.g.vcf.gz')
    samples_df.to_csv(samples_path, index=False, sep='\t', na_rep='NA')
    logger.info(f'Saved metadata with updated GVCFs to {samples_path}')
    combined_mt_path = join(combiner_bucket, 'genomes.mt')

    # if not utils.file_exists(combined_mt_path):
    #     combiner_job = dataproc.hail_dataproc_job(
    #         b,
    #         f'scripts/run_python_script.py '
    #         f'combine_gvcfs.py '
    #         f'--meta-csv {samples_path} '
    #         f'--out-mt {combined_mt_path} '
    #         f'--bucket {combiner_bucket}/work '
    #         f'--hail-billing {billing_project} ',
    #         max_age='8h',
    #         packages=utils.DATAPROC_PACKAGES,
    #         num_secondary_workers=10,
    #         depends_on=subset_gvcf_jobs,
    #         job_name='Combine GVCFs',
    #     )
    # else:
    #     combiner_job = b.new_job('Combine GVCFs')

    hard_filtered_samples_ht_path = join(combiner_bucket, 'hard_filters.ht')
    meta_ht_path = join(combiner_bucket, 'meta.ht')
    # if not any(
    #     utils.file_exists(fp) for fp in [hard_filtered_samples_ht_path, meta_ht_path]
    # ):
    #     age_csv = join(base_bucket, 'age.csv')
    #     if utils.file_exists(age_csv):
    #         age_csv_param = f'--age-csv {age_csv} '
    #     else:
    #         age_csv_param = ''
    #     sample_qc_job = dataproc.hail_dataproc_job(
    #         b,
    #         f'scripts/run_python_script.py '
    #         f'sample_qc.py --overwrite '
    #         f'--mt {combined_mt_path} '
    #         f'{age_csv_param}'
    #         f'--meta-csv {samples_path} '
    #         f'--bucket {combiner_bucket} '
    #         f'--out-hardfiltered-samples-ht {hard_filtered_samples_ht_path} '
    #         f'--out-meta-ht {meta_ht_path} '
    #         f'--hail-billing {billing_project} ',
    #         max_age='8h',
    #         packages=utils.DATAPROC_PACKAGES,
    #         num_secondary_workers=10,
    #         depends_on=[combiner_job],
    #         job_name='Sample QC',
    #     )
    # else:
    #     sample_qc_job = b.new_job('Sample QC')

    variant_qc_bucket = join(output_bucket, 'variant_qc')
    freq_ht_path = join(variant_qc_bucket, 'frequencies.ht')
    info_ht_path = join(variant_qc_bucket, 'info.ht')
    allele_data_ht_path = join(variant_qc_bucket, 'allele_data.ht')
    qc_ac_ht_path = join(variant_qc_bucket, 'qc_ac.ht')
    # vep_ht = join(variant_qc_bucket, 'vep.ht')
    rf_result_ht_path = join(variant_qc_bucket, 'rf_result.ht')
    # if not all(
    #     utils.file_exists(fp)
    #     for fp in [info_ht_path, allele_data_ht_path, qc_ac_ht_path]
    # ):
    #     rf_anno_job = dataproc.hail_dataproc_job(
    #         b,
    #         f'scripts/run_python_script.py '
    #         f'generate_qc_annotations.py --reuse '
    #         f'--split-multiallelic '
    #         f'--mt {combined_mt_path} '
    #         f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
    #         f'--meta-ht {meta_ht_path} '
    #         f'--out-info-ht {info_ht_path} '
    #         f'--out-allele-data-ht {allele_data_ht_path} '
    #         f'--out-qc-ac-ht {qc_ac_ht_path} '
    #         f'--bucket {variant_qc_bucket} ',
    #         max_age='8h',
    #         packages=utils.DATAPROC_PACKAGES,
    #         num_secondary_workers=10,
    #         depends_on=[sample_qc_job],
    #         job_name='RF: gen QC anno',
    #         vep='GRCh38',
    #     )
    # else:
    #     rf_anno_job = b.new_job('RF: gen QC anno')

    # if not utils.file_exists(freq_ht_path):
    #     rf_freq_data_job = dataproc.hail_dataproc_job(
    #         b,
    #         f'scripts/run_python_script.py '
    #         f'generate_freq_data.py --reuse '
    #         f'--mt {combined_mt_path} '
    #         f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
    #         f'--meta-ht {meta_ht_path} '
    #         f'--out-ht {freq_ht_path} '
    #         f'--bucket {variant_qc_bucket} ',
    #         max_age='8h',
    #         packages=utils.DATAPROC_PACKAGES,
    #         num_secondary_workers=10,
    #         depends_on=[sample_qc_job],
    #         job_name='RF: gen freq data',
    #     )
    # else:
    #     rf_freq_data_job = b.new_job('RF: gen freq data')
    
    if not utils.file_exists(rf_result_ht_path):
        rf_job = dataproc.hail_dataproc_job(
            b,
            f'scripts/run_python_script.py '
            f'random_forest.py --reuse '
            f'--info-ht {info_ht_path} '
            f'--freq-ht {freq_ht_path} '
            f'--allele-data-ht {allele_data_ht_path} '
            f'--qc-ac-ht {qc_ac_ht_path} '
            f'--bucket {variant_qc_bucket} '
            '--use-adj-genotypes '
            f'--out-ht {rf_result_ht_path} ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            # depends_on=[rf_freq_data_job, rf_anno_job],
            job_name='RF: run',
        )
    else:
        rf_job = b.new_job('RF: run')

    rf_job.always_run()

    # make_vqsr_jobs(
    #     b,
    #     combined_mt_path=combined_mt_path,
    #     gvcf_count=len(gvcfs),
    #     work_bucket=join(output_bucket, 'vqsr'),
    #     callset_name=callset_name,
    #     is_small_callset=is_small_callset,
    #     is_huge_callset=is_huge_callset,
    #     depends_on=[combiner_job],
    #     skip_allele_specific_annotations=skip_allele_specific_annotations,
    #     snp_filter_level=snp_filter_level,
    #     indel_filter_level=indel_filter_level,
    # )

    b.run(dry_run=dry_run, delete_scratch_on_exit=not keep_scratch)


def add_reblock_gvcfs_step(
    b: hb.Batch,
    input_gvcf: hb.ResourceGroup,
) -> Job:
    """
    Runs ReblockGVCFs to annotate with allele-specific VCF INFO fields
    required for recalibration
    """
    j = b.new_job('ReblocGVCFs')
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
        --drop-low-quals \\
        -do-qual-approx \\
        -O {j.output_gvcf['g.vcf.gz']} \\
        --create-output-variant-index true"""
    )
    return j


def add_subset_noalt_step(
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
    if utils.file_exists(output_gvcf_path):
        return j

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
