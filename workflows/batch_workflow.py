"""
Hail Batch workflow to perform VQSR on a WGS germline callset.
Compilation from the following two WDL workflows:
1. hail-ukbb-200k-callset/GenotypeAndFilter.AS.wdl
2. The Broad VQSR workflow:
   https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/dna_seq/germline/joint_genotyping/JointGenotyping.wdl
   documented here:
   https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
Translated from WDL with a help of Janis:
https://github.com/PMCC-BioinformaticsCore/janis

The input is a VCF file, build the following way:
1. called variants in single samples with GATK HaplotypeCaller
2. post-process GVCFs with GATK ReblockGVCF
3. combine the GVCFs with Hail combiner into a Hail Matrix Table
   using `scripts/combine_gvcfs.py`
4. optionally, perform sample-level filtering with `scripts/sample_qc.py`
5. export the Matrix Table to a multi-sample site-only VCF with
   `scripts/mt_to_vcf.py`

The output is a VCF file <output_bucket>/<callset>-recalibrated.vcf.gz,
as well as a QC file <output_bucket>/<callset>-eval.txt
and R scripts to plot VQSR models: <output_bucket>/plot-snps-recal.Rscript
and <output_bucket>/plot-indels-recal.Rscript
"""

import os
from os.path import join
from typing import List
import logging
import click
import hailtop.batch as hb
from hailtop.batch.job import Job
from analysis_runner import dataproc

from joint_calling import utils

logger = logging.getLogger('joint-calling')
logger.setLevel('INFO')


GATK_VERSION = '4.2.0.0'
GATK_DOCKER = f'us.gcr.io/broad-gatk/gatk:{GATK_VERSION}'
# GnarlyGenotyper crashes with NullPointerException when using GATK docker
GNARLY_DOCKER = 'gcr.io/broad-dsde-methods/gnarly_genotyper:hail_ukbb_300K'
DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver'

BROAD_REF_BUCKET = 'gs://gcp-public-data--broad-references/hg38/v0'

DATAPROC_PACKAGES = [
    'joint-calling',
    'click',
    'cpg-gnomad',
    'google',
    'slackclient',
    'fsspec',
    'sklearn',
    'gcloud',
]


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
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option(
    '--billing-project',
    'billing_project',
    type=str,
    default=os.getenv('HAIL_BILLING_PROJECT'),
)
# Reference files
@click.option(
    '--unpadded_intervals_file',
    'unpadded_intervals_file',
    type=str,
    default=os.path.join(BROAD_REF_BUCKET, 'hg38.even.handcurated.20k.intervals'),
)
@click.option(
    '--ref_fasta',
    'ref_fasta',
    type=str,
    default=os.path.join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.fasta'),
)
@click.option(
    '--ref_fasta_index',
    'ref_fasta_index',
    type=str,
    default=os.path.join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.fasta.fai'),
)
@click.option(
    '--ref_dict',
    'ref_dict',
    type=str,
    default=os.path.join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.dict'),
)
@click.option(
    '--dbsnp_vcf',
    'dbsnp_vcf',
    type=str,
    default=os.path.join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf'),
)
@click.option(
    '--dbsnp_vcf_index',
    'dbsnp_vcf_index',
    type=str,
    default=os.path.join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf.idx'),
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
    '--eval_interval_list',
    'eval_interval_list',
    type=str,
    default=os.path.join(BROAD_REF_BUCKET, 'wgs_evaluation_regions.hg38.interval_list'),
)
@click.option(
    '--hapmap_resource_vcf',
    'hapmap_resource_vcf',
    type=str,
    default=os.path.join(BROAD_REF_BUCKET, 'hapmap_3.3.hg38.vcf.gz'),
)
@click.option(
    '--hapmap_resource_vcf_index',
    'hapmap_resource_vcf_index',
    type=str,
    default=os.path.join(BROAD_REF_BUCKET, 'hapmap_3.3.hg38.vcf.gz.tbi'),
)
@click.option(
    '--omni_resource_vcf',
    'omni_resource_vcf',
    type=str,
    default=os.path.join(BROAD_REF_BUCKET, '1000G_omni2.5.hg38.vcf.gz'),
)
@click.option(
    '--omni_resource_vcf_index',
    'omni_resource_vcf_index',
    type=str,
    default=os.path.join(BROAD_REF_BUCKET, '1000G_omni2.5.hg38.vcf.gz.tbi'),
)
@click.option(
    '--one_thousand_genomes_resource_vcf',
    'one_thousand_genomes_resource_vcf',
    type=str,
    default=os.path.join(
        BROAD_REF_BUCKET, '1000G_phase1.snps.high_confidence.hg38.vcf.gz'
    ),
)
@click.option(
    '--one_thousand_genomes_resource_vcf_index',
    'one_thousand_genomes_resource_vcf_index',
    type=str,
    default=os.path.join(
        BROAD_REF_BUCKET, '1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi'
    ),
)
@click.option(
    '--mills_resource_vcf',
    'mills_resource_vcf',
    type=str,
    default=os.path.join(
        BROAD_REF_BUCKET, 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
    ),
)
@click.option(
    '--mills_resource_vcf_index',
    'mills_resource_vcf_index',
    type=str,
    default=os.path.join(
        BROAD_REF_BUCKET, 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'
    ),
)
@click.option(
    '--axiom_poly_resource_vcf',
    'axiom_poly_resource_vcf',
    type=str,
    default=os.path.join(
        BROAD_REF_BUCKET, 'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz'
    ),
)
@click.option(
    '--axiom_poly_resource_vcf_index',
    'axiom_poly_resource_vcf_index',
    type=str,
    default=os.path.join(
        BROAD_REF_BUCKET,
        'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi',
    ),
)
@click.option('--dbsnp_resource_vcf', 'dbsnp_resource_vcf', type=str)
@click.option('--dbsnp_resource_vcf_index', 'dbsnp_resource_vcf_index', type=str)
# ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
# than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
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
    dry_run: bool,
    billing_project: str,
    unpadded_intervals_file: str,
    ref_fasta: str,
    ref_fasta_index: str,
    ref_dict: str,
    dbsnp_vcf: str,
    dbsnp_vcf_index: str,
    snp_recalibration_tranche_values: List[float],
    snp_recalibration_annotation_values: List[str],
    indel_recalibration_tranche_values: List[float],
    indel_recalibration_annotation_values: List[str],
    eval_interval_list: str,
    hapmap_resource_vcf: str,
    hapmap_resource_vcf_index: str,
    omni_resource_vcf: str,
    omni_resource_vcf_index: str,
    one_thousand_genomes_resource_vcf: str,
    one_thousand_genomes_resource_vcf_index: str,
    mills_resource_vcf: str,
    mills_resource_vcf_index: str,
    axiom_poly_resource_vcf: str,
    axiom_poly_resource_vcf_index: str,
    dbsnp_resource_vcf: str,
    dbsnp_resource_vcf_index: str,
    excess_het_threshold: float,
    snp_filter_level: float,
    indel_filter_level: float,
    skip_allele_specific_annotations: bool,
):
    """
    Run a Hail Batch workflow that performs the VQSR filtering on a WGS
    germline callset
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
    output_bucket = f'gs://cpg-{callset_name}-{output_bucket_suffix}/joint_vcf/{callset_version}/work'

    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=join(output_bucket, 'hail').replace('gs://', ''),
    )
    b = hb.Batch('Joint Calling', backend=backend)

    os.system('echo "HAIL_QUERY_BACKEND=$HAIL_QUERY_BACKEND"')
    os.environ['HAIL_QUERY_BACKEND'] = 'spark'
    os.system('echo "HAIL_QUERY_BACKEND=$HAIL_QUERY_BACKEND"')
    assert os.environ.get('HAIL_QUERY_BACKEND') != 'service'
    utils.init_hail('joint-calling-workflow')

    # TODO: merge with existing data
    # TODO: fix impute_type

    samples_df = utils.find_inputs(input_buckets, skip_qc=skip_qc)
    samples_path = join(output_bucket, 'samples.csv')
    if not utils.file_exists(samples_path):
        samples_df.to_csv(samples_path, index=False, sep='\t', na_rep='NA')
    logger.info(f'Saved metadata to {samples_path}')

    # logger.info(f'Testing reading...')
    # df = pd.read_table(samples_path)
    # print(df)
    # input_meta_ht = hl.Table.from_pandas(df)
    # input_meta_ht.show()

    gvcfs = [
        b.read_input_group(**{'g.vcf.gz': gvcf, 'g.vcf.gz.tbi': gvcf + '.tbi'})
        for gvcf in list(samples_df.gvcf)
    ]

    # Make a 2.5:1 interval number to samples in callset ratio interval list.
    # We allow overriding the behavior by specifying the desired number of vcfs
    # to scatter over for testing / special requests.
    scatter_count_scale_factor = 0.15
    scatter_count = int(round(scatter_count_scale_factor * len(gvcfs)))
    scatter_count = max(scatter_count, 2)

    is_small_callset = len(gvcfs) < 1000
    # 1. For small callsets, we don't apply the ExcessHet filtering.
    # 2. For small callsets, we gather the VCF shards and collect QC metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics
    # collected from them.
    is_huge_callset = len(gvcfs) >= 100000
    # For huge callsets, we allocate more memory for the SNPs Create Model step

    small_disk = 30 if is_small_callset else (50 if not is_huge_callset else 100)
    medium_disk = 50 if is_small_callset else (100 if not is_huge_callset else 200)
    huge_disk = 100 if is_small_callset else (500 if not is_huge_callset else 2000)

    ref_fasta = b.read_input_group(
        base=ref_fasta,
        dict=ref_dict
        or (
            ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
            + '.dict'
        ),
        fai=ref_fasta_index or (ref_fasta + '.fai'),
    )
    dbsnp_vcf = b.read_input_group(base=dbsnp_vcf, index=dbsnp_vcf_index)
    eval_interval_list = b.read_input(eval_interval_list)
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
    dbsnp_resource_vcf = (
        b.read_input_group(base=dbsnp_resource_vcf, index=dbsnp_resource_vcf_index)
        if dbsnp_resource_vcf
        else dbsnp_vcf
    )
    # pylint: disable=unused-variable
    noalt_regions = b.read_input('gs://cpg-reference/hg38/v0/noalt.bed')

    # reblocked_gvcfs = [
    #     add_reblock_gvcfs_step(b, gvcf, small_disk).output_gvcf for gvcf in gvcfs
    # ]
    combiner_bucket = os.path.join(output_bucket, 'combiner')
    # combiner_gvcf_bucket = os.path.join(output_bucket, 'combiner', 'gvcfs')
    # subset_gvcf_jobs = [
    #     add_subset_noalt_step(
    #         b,
    #         input_gvcf=gvcf,
    #         output_gvcf_path=join(combiner_gvcf_bucket, sample + '.g.vcf.gz'),
    #         disk_size=small_disk,
    #         noalt_regions=noalt_regions,
    #     )
    #     for sample, gvcf in zip(list(samples_df.s), reblocked_gvcfs)
    # ]

    combined_mt_path = join(combiner_bucket, 'genomes.mt')
    hard_filtered_samples_ht_path = join(combiner_bucket, 'hard_filters.ht')
    meta_ht_path = join(combiner_bucket, 'meta.ht')
    combined_vcf_path = join(combiner_bucket, 'combined.vcf.gz')
    vcf_buckets_cmdl = ' '.join([f'--bucket-with-vcfs {ib}' for ib in input_buckets])
    combiner_job = dataproc.hail_dataproc_job(
        b,
        f'run_python_script.py '
        f'combine_gvcfs.py --reuse '
        f'{vcf_buckets_cmdl} '
        f'{"--skip-qc " if skip_qc else ""}'
        f'--out-mt {combined_mt_path} '
        f'--bucket {combiner_bucket}/work '
        f'--hail-billing {billing_project} ',
        max_age='8h',
        packages=DATAPROC_PACKAGES,
        num_secondary_workers=10,
        # depends_on=subset_gvcf_jobs,
        job_name='Combine GVCFs',
    )
    sample_qc_job = dataproc.hail_dataproc_job(
        b,
        f'run_python_script.py '
        f'sample_qc.py '
        f'--mt {combined_mt_path} '
        f'--meta-csv {samples_path} '
        f'--bucket {combiner_bucket} '
        f'--out-hardfiltered-samples-ht {hard_filtered_samples_ht_path} '
        f'--out-meta-ht {meta_ht_path} '
        f'--hail-billing {billing_project} ',
        max_age='8h',
        packages=DATAPROC_PACKAGES,
        num_secondary_workers=10,
        depends_on=[combiner_job],
        job_name='Sample QC',
    )
    mt_to_vcf_job = dataproc.hail_dataproc_job(
        b,
        f'run_python_script.py '
        f'mt_to_vcf.py --overwrite '
        f'--mt {combined_mt_path} '
        f'-o {combined_vcf_path} ',
        max_age='8h',
        packages=DATAPROC_PACKAGES,
        num_secondary_workers=10,
        depends_on=[combiner_job],
        job_name='MT to VCF',
    )

    variant_qc_bucket = join(output_bucket, 'variant_qc')
    freq_ht_path = join(variant_qc_bucket, 'frequencies.ht')
    info_ht_path = join(variant_qc_bucket, 'info.ht')
    allele_data_ht_path = join(variant_qc_bucket, 'allele_data.ht')
    qc_ac_ht_path = join(variant_qc_bucket, 'qc_ac.ht')
    # vep_ht = join(variant_qc_bucket, 'vep.ht')
    rf_result_ht_path = join(variant_qc_bucket, 'rf_result.ht')
    rf_anno_job = dataproc.hail_dataproc_job(
        b,
        f'run_python_script.py '
        f'generate_qc_annotations.py --split-multiallelic --overwrite '
        f'--mt {combined_mt_path} '
        f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
        f'--meta-ht {meta_ht_path} '
        f'--out-info-ht {info_ht_path} '
        f'--out-allele-data-ht {allele_data_ht_path} '
        f'--out-qc-ac-ht {qc_ac_ht_path} '
        f'--bucket {combiner_bucket} ',
        max_age='8h',
        packages=DATAPROC_PACKAGES,
        num_secondary_workers=10,
        depends_on=[sample_qc_job],
        job_name='RF: gen QC anno',
    )
    rf_freq_data_job = dataproc.hail_dataproc_job(
        b,
        f'run_python_script.py '
        f'generate_freq_data.py --overwrite '
        f'--mt {combined_mt_path} '
        f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
        f'--meta-ht {meta_ht_path} '
        f'--out-ht {freq_ht_path} '
        f'--bucket {combiner_bucket} ',
        max_age='8h',
        packages=DATAPROC_PACKAGES,
        num_secondary_workers=10,
        depends_on=[rf_anno_job],
        job_name='RF: gen freq data',
    )
    rf_job = dataproc.hail_dataproc_job(
        b,
        f'run_python_script.py '
        f'random_forest.py --overwrite '
        f'--info-ht {info_ht_path} '
        f'--freq-ht {freq_ht_path} '
        f'--allele-data-ht {allele_data_ht_path} '
        f'--qc-ac-ht {qc_ac_ht_path} '
        f'--out-ht {rf_result_ht_path} '
        f'--bucket {combiner_bucket} '
        f'-o {rf_result_ht_path} ',
        max_age='8h',
        packages=DATAPROC_PACKAGES,
        num_secondary_workers=10,
        depends_on=[rf_freq_data_job],
        job_name='RF: main',
    )
    rf_job.always_run()

    split_intervals_job = add_split_intervals_step(
        b,
        unpadded_intervals_file,
        scatter_count,
        ref_fasta,
        disk_size=small_disk,
    )
    intervals = split_intervals_job.intervals

    tabix_job = add_tabix_step(b, combined_vcf_path, medium_disk)
    tabix_job.depends_on(mt_to_vcf_job)

    gnarly_output_vcfs = [
        add_gnarly_genotyper_on_vcf_step(
            b,
            combined_gvcf=tabix_job.combined_vcf,
            interval=intervals[f'interval_{idx}'],
            ref_fasta=ref_fasta,
            dbsnp_vcf=dbsnp_vcf,
            disk_size=medium_disk,
        ).output_vcf
        for idx in range(scatter_count)
    ]

    if not is_small_callset:
        # ExcessHet filtering applies only to callsets with a large number of samples,
        # e.g. hundreds of unrelated samples. Small cohorts should not trigger ExcessHet
        # filtering as values should remain small. Note cohorts of consanguinous samples
        # will inflate ExcessHet, and it is possible to limit the annotation to founders
        # for such cohorts by providing a pedigree file during variant calling.
        hard_filtered_vcfs = [
            add_hard_filter_step(
                b,
                input_vcf=gnarly_output_vcfs[idx],
                excess_het_threshold=excess_het_threshold,
                disk_size=medium_disk,
            ).output_vcf
            for idx in range(scatter_count)
        ]
    else:
        hard_filtered_vcfs = gnarly_output_vcfs
    # hard_filtered_vcfs = [
    #     b.read_input_group(**{
    #         'vcf.gz': f'gs://playground-au/batch/859e9a/{idx + 2}/output_vcf.vcf.gz',
    #         'vcf.gz.tbi': f'gs://playground-au/batch/859e9a/{idx + 2}/output_vcf.vcf.gz.tbi'})
    #     for idx in range(scatter_count)
    # ]

    sites_only_vcfs = [
        add_make_sites_only_vcf_step(
            b,
            input_vcf=hard_filtered_vcfs[idx],
            disk_size=medium_disk,
        ).sites_only_vcf
        for idx in range(scatter_count)
    ]
    # sites_only_vcfs = [
    #     b.read_input_group(**{
    #         'vcf.gz': f'gs://playground-au/batch/859e9a/{idx + 9}/'
    #             'sites_only_vcf.vcf.gz',
    #         'vcf.gz.tbi': f'gs://playground-au/batch/859e9a/{idx + 9}/'
    #             'sites_only_vcf.vcf.gz.tbi'})
    #     for idx in range(scatter_count)
    # ]

    sites_only_gathered_vcf = add_sites_only_gather_vcf_step(
        b,
        input_vcfs=sites_only_vcfs,
        disk_size=medium_disk,
    ).output_vcf
    # sites_only_gathered_vcf = b.read_input_group(**{
    #     'vcf.gz': 'gs://playground-au/batch/1e0bc3/1/output_vcf.vcf.gz',
    #     'vcf.gz.tbi': 'gs://playground-au/batch/1e0bc3/1/output_vcf.vcf.gz.tbi'
    # })
    #
    indels_variant_recalibrator_job = add_indels_variant_recalibrator_step(
        b,
        sites_only_variant_filtered_vcf=sites_only_gathered_vcf,
        recalibration_tranche_values=indel_recalibration_tranche_values,
        recalibration_annotation_values=indel_recalibration_annotation_values,
        mills_resource_vcf=mills_resource_vcf,
        axiom_poly_resource_vcf=axiom_poly_resource_vcf,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        use_allele_specific_annotations=not skip_allele_specific_annotations,
        disk_size=small_disk,
        output_bucket=output_bucket,
    )
    indels_recalibration = indels_variant_recalibrator_job.recalibration
    indels_tranches = indels_variant_recalibrator_job.tranches
    # indels_recalibration = b.read_input_group(
    #     base='gs://playground-au/batch/859e9a/17/recalibration',
    #     index='gs://playground-au/batch/859e9a/17/recalibration.idx',
    # )
    # indels_tranches = b.read_input('gs://playground-au/batch/859e9a/17/tranches')

    snp_max_gaussians = 6
    if is_small_callset:
        snp_max_gaussians = 4
    elif is_huge_callset:
        snp_max_gaussians = 8

    if is_huge_callset:
        # Run SNP recalibrator in a scattered mode
        model_file = add_snps_variant_recalibrator_create_model_step(
            b,
            sites_only_variant_filtered_vcf=sites_only_gathered_vcf,
            recalibration_tranche_values=snp_recalibration_tranche_values,
            recalibration_annotation_values=snp_recalibration_annotation_values,
            hapmap_resource_vcf=hapmap_resource_vcf,
            omni_resource_vcf=omni_resource_vcf,
            one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
            dbsnp_resource_vcf=dbsnp_resource_vcf,
            disk_size=small_disk,
            output_bucket=output_bucket,
            use_allele_specific_annotations=not skip_allele_specific_annotations,
            is_small_callset=is_small_callset,
            is_huge_callset=is_huge_callset,
            max_gaussians=snp_max_gaussians,
        ).model_file
        # model_file = b.read_input('gs://playground-au/batch/859e9a/18/model_report')

        snps_recalibrator_jobs = [
            add_snps_variant_recalibrator_scattered_step(
                b,
                sites_only_variant_filtered_vcf=sites_only_vcfs[idx],
                recalibration_tranche_values=snp_recalibration_tranche_values,
                recalibration_annotation_values=snp_recalibration_annotation_values,
                model_file=model_file,
                hapmap_resource_vcf=hapmap_resource_vcf,
                omni_resource_vcf=omni_resource_vcf,
                one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
                dbsnp_resource_vcf=dbsnp_resource_vcf,
                disk_size=small_disk,
                max_gaussians=snp_max_gaussians,
                use_allele_specific_annotations=not skip_allele_specific_annotations,
            )
            for idx in range(len(sites_only_vcfs))
        ]
        snps_recalibrations = [j.recalibration for j in snps_recalibrator_jobs]
        snps_tranches = [j.tranches for j in snps_recalibrator_jobs]
        # snp_tranches = [
        #     b.read_input(f'gs://playground-au/batch/df311d/{idx + 1}/tranches')
        #     for idx in range(scatter_count)
        # ]
        # snp_recalibrations = [
        #     b.read_input(f'gs://playground-au/batch/df311d/{idx + 1}/recalibration')
        #     for idx in range(scatter_count)
        # ]
        snps_gathered_tranches = add_snps_gather_tranches_step(
            b,
            tranches=snps_tranches,
            disk_size=small_disk,
        ).out_tranches

        recalibrated_vcfs = [
            add_apply_recalibration_step(
                b,
                input_vcf=hard_filtered_vcfs[idx],
                indels_recalibration=indels_recalibration,
                indels_tranches=indels_tranches,
                snps_recalibration=snps_recalibrations[idx],
                snps_tranches=snps_gathered_tranches,
                disk_size=medium_disk,
                use_allele_specific_annotations=not skip_allele_specific_annotations,
                indel_filter_level=indel_filter_level,
                snp_filter_level=snp_filter_level,
            ).recalibrated_vcf
            for idx in range(len(hard_filtered_vcfs))
        ]
        final_gathered_vcf = add_final_gather_vcf_step(
            b,
            input_vcfs=recalibrated_vcfs,
            disk_size=huge_disk,
            output_vcf_path=os.path.join(
                output_bucket, callset_name + '-recalibrated.vcf.gz'
            ),
        ).output_vcf

    else:
        snps_recalibrator_job = add_snps_variant_recalibrator_step(
            b,
            sites_only_variant_filtered_vcf=sites_only_gathered_vcf,
            recalibration_tranche_values=snp_recalibration_tranche_values,
            recalibration_annotation_values=snp_recalibration_annotation_values,
            hapmap_resource_vcf=hapmap_resource_vcf,
            omni_resource_vcf=omni_resource_vcf,
            one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
            dbsnp_resource_vcf=dbsnp_resource_vcf,
            disk_size=small_disk,
            max_gaussians=snp_max_gaussians,
            use_allele_specific_annotations=not skip_allele_specific_annotations,
            output_bucket=output_bucket,
        )
        snps_recalibration = snps_recalibrator_job.recalibration
        snps_tranches = snps_recalibrator_job.tranches

        gathered_vcf = add_final_gather_vcf_step(
            b,
            input_vcfs=hard_filtered_vcfs,
            disk_size=huge_disk,
        ).output_vcf

        final_gathered_vcf = add_apply_recalibration_step(
            b,
            input_vcf=gathered_vcf,
            indels_recalibration=indels_recalibration,
            indels_tranches=indels_tranches,
            snps_recalibration=snps_recalibration,
            snps_tranches=snps_tranches,
            disk_size=medium_disk,
            use_allele_specific_annotations=not skip_allele_specific_annotations,
            indel_filter_level=indel_filter_level,
            snp_filter_level=snp_filter_level,
        ).recalibrated_vcf

    add_variant_eval_step(
        b,
        input_vcf=final_gathered_vcf,
        ref_fasta=ref_fasta,
        dbsnp_vcf=dbsnp_vcf,
        output_path=os.path.join(output_bucket, callset_name + '-eval.txt'),
        disk_size=huge_disk,
    )

    b.run(dry_run=dry_run, delete_scratch_on_exit=not keep_scratch)


def add_reblock_gvcfs_step(
    b: hb.Batch,
    input_gvcf: hb.ResourceGroup,
    disk_size: int,
) -> Job:
    """
    Runs ReblockGVCFs to annotate with allele-specific VCF INFO fields
    required for recalibration
    """
    j = b.new_job('ReblocGVCFs')
    j.image(GATK_DOCKER)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'{disk_size}G')
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
    disk_size: int,
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
    j.storage(f'{disk_size}G')
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


def add_tabix_step(
    b: hb.Batch,
    vcf_path: str,
    disk_size: int,
) -> Job:
    """
    Regzip and tabix the combined VCF (for some reason the one output with mt2vcf
    is not block-gzipped)
    """
    j = b.new_job('Tabix')
    j.image('quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2')
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
    j = b.new_job('SplitIntervals')
    j.image(GATK_DOCKER)
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

    gatk --java-options -Xms{mem_gb - 1}g SplitIntervals \\
      -L {interval_list} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {ref_fasta.base} \\
      -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
      """
    )
    return j


def add_gnarly_genotyper_on_vcf_step(
    b: hb.Batch,
    combined_gvcf: hb.ResourceGroup,
    interval: hb.ResourceGroup,
    ref_fasta: hb.ResourceGroup,
    dbsnp_vcf: hb.ResourceGroup,
    disk_size: int,
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
    j = b.new_job('GnarlyGenotyperOnVcf')
    # GnarlyGenotyper crashes with NullPointerException when using standard GATK docker
    j.image(GNARLY_DOCKER)
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
      -L {interval} \\
      --create-output-variant-index"""
    )
    return j


def add_hard_filter_step(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    excess_het_threshold: float,
    disk_size: int,
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
    j = b.new_job('HardFilter')
    j.image(GATK_DOCKER)
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
) -> Job:
    """
    Create sites-only VCF with only site-level annotations.
    Speeds up the analysis in the modeling step.

    Returns: a Job object with a single output j.sites_only_vcf of type ResourceGroup
    """
    j = b.new_job('MakeSitesOnlyVcf')
    j.image(GATK_DOCKER)
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
      -O {j.sites_only_vcf['vcf.gz']}"""
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
    j = b.new_job('SitesOnlyGatherVcf')
    j.image(GATK_DOCKER)
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
    recalibration_tranche_values: List[float],
    recalibration_annotation_values: List[str],
    mills_resource_vcf: hb.ResourceGroup,
    axiom_poly_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    use_allele_specific_annotations: bool,
    disk_size: int,
    output_bucket: str = None,
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
    j = b.new_job('IndelsVariantRecalibrator')
    j.image(GATK_DOCKER)
    mem_gb = 32
    j.memory(f'{mem_gb}G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in recalibration_tranche_values])
    an_cmdl = ' '.join([f'-an {v}' for v in recalibration_annotation_values])
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
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''} \\
      --max-gaussians {max_gaussians} \\
      -resource:mills,known=false,training=true,truth=true,prior=12 {mills_resource_vcf.base} \\
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {axiom_poly_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 {dbsnp_resource_vcf.base} \\
      --rscript-file {j.indel_rscript_file}"""
    )
    if output_bucket:
        b.write_output(
            j.indel_rscript_file,
            os.path.join(output_bucket, 'plot-indels-recal.Rscript'),
        )
    return j


def add_snps_variant_recalibrator_create_model_step(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    recalibration_tranche_values: List[float],
    recalibration_annotation_values: List[str],
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    output_bucket: str = None,
    use_allele_specific_annotations: bool = True,
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
    The latter is usedful to produce the optional tranche plot.
    """
    j = b.new_job('SNPsVariantRecalibratorCreateModel')
    j.image(GATK_DOCKER)
    mem_gb = 64 if not is_small_callset else 128
    j.memory(f'{mem_gb}G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    downsample_factor = 75 if is_huge_callset else 10

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in recalibration_tranche_values])
    an_cmdl = ' '.join([f'-an {v}' for v in recalibration_annotation_values])
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
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''} \\
      --sample-every-Nth-variant {downsample_factor} \\
      --output-model {j.model_file} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base} \\
      --rscript-file {j.snp_rscript}

      ln {j.snp_rscript}.pdf {j.snp_rscript_pdf}"""
    )
    if output_bucket:
        b.write_output(
            j.snp_rscript_pdf,
            os.path.join(output_bucket, 'recalibration-indels-features.pdf'),
        )
    return j


def add_snps_variant_recalibrator_scattered_step(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    model_file: hb.ResourceGroup,
    recalibration_tranche_values: List[float],
    recalibration_annotation_values: List[str],
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    use_allele_specific_annotations: bool = True,
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
    j = b.new_job('SNPsVariantRecalibratorScattered')

    j.image(GATK_DOCKER)
    mem_gb = 64  # ~ twice the sum of all input resources and input VCF sizes
    j.memory(f'{mem_gb}G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in recalibration_tranche_values])
    an_cmdl = ' '.join([f'-an {v}' for v in recalibration_annotation_values])
    j.command(
        f"""set -euo pipefail

    MODEL_REPORT={model_file}

    gatk --java-options -Xms{mem_gb - 1}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''} \\
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
    recalibration_tranche_values: List[float],
    recalibration_annotation_values: List[str],
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    output_bucket: str,
    disk_size: int,
    use_allele_specific_annotations: bool = True,
    max_gaussians: int = 4,
) -> Job:
    """
    Recalibrate SNPs in one run (alternative to scatter-gather approach)
    """
    j = b.new_job('SNPsVariantRecalibrator')

    j.image(GATK_DOCKER)
    mem_gb = 64  # ~ twice the sum of all input resources and input VCF sizes
    j.memory(f'{mem_gb}G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in recalibration_tranche_values])
    an_cmdl = ' '.join([f'-an {v}' for v in recalibration_annotation_values])
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
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base} \\
      --rscript-file {j.snp_rscript}

      ln {j.snp_rscript}.pdf {j.snp_rscript_pdf}
      ln {j.tranches}.pdf {j.tranches_pdf}"""
    )

    if output_bucket:
        b.write_output(
            j.snp_rscript_pdf,
            os.path.join(output_bucket, 'recalibration-snps-features.pdf'),
        )
        b.write_output(
            j.tranches_pdf,
            os.path.join(output_bucket, 'recalibration-snps-tranches.pdf'),
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
    j = b.new_job('SNPGatherTranches')
    j.image(GATK_DOCKER)
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
    indel_filter_level: float = 99.0,
    snp_filter_level: float = 99.7,
    use_allele_specific_annotations: bool = True,
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
    j = b.new_job('ApplyRecalibration')
    j.image(GATK_DOCKER)
    j.memory('8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        recalibrated_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O tmp.indel.recalibrated.vcf \\
      -V {input_vcf['vcf.gz']} \\
      --recal-file {indels_recalibration} \\
      --tranches-file {indels_tranches} \\
      --truth-sensitivity-filter-level {indel_filter_level} \\
      --create-output-variant-index true \\
      -mode INDEL \\
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''}

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O {j.recalibrated_vcf['vcf.gz']} \\
      -V tmp.indel.recalibrated.vcf \\
      --recal-file {snps_recalibration} \\
      --tranches-file {snps_tranches} \\
      --truth-sensitivity-filter-level {snp_filter_level} \\
      --create-output-variant-index true \\
      -mode SNP \\
      {'--use-allele-specific-annotations' if use_allele_specific_annotations else ''}"""
    )
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
    j = b.new_job('CollectMetricsSharded')
    j.image(GATK_DOCKER)
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


def add_final_gather_vcf_step(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceGroup],
    disk_size: int,
    output_vcf_path: str = None,
) -> Job:
    """
    Combines recalibrated VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    j = b.new_job('FinalGatherVcf')
    j.image(GATK_DOCKER)
    j.memory(f'8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
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
      --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}"""
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j


def add_variant_eval_step(
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
    j = b.new_job('VariantEval')
    j.image(GATK_DOCKER)
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
    j = b.new_job('GatherVariantCallingMetrics')
    j.image(GATK_DOCKER)
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


if __name__ == '__main__':
    main()  # pylint: disable=E1120
