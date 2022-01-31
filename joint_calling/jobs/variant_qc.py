"""
Variant QC Hail-query jobs
"""

from os.path import join
from typing import List, Optional, Dict
import logging
import hailtop.batch as hb
from hailtop.batch.job import Job

from joint_calling import utils
from joint_calling.dataproc import add_job
from joint_calling.jobs.vqsr import add_vqsr_jobs

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def add_variant_qc_jobs(
    b: hb.Batch,
    work_bucket: str,
    web_bucket: str,
    raw_combined_mt_path: str,
    hard_filter_ht_path: str,
    meta_ht_path: str,
    out_filtered_combined_mt_path: str,
    out_filtered_vcf_ptrn_path: str,
    sample_count: int,
    ped_file: Optional[str],
    overwrite: bool,
    vqsr_params_d: Dict,
    scatter_count: int,
    is_test: bool,
    project_name: str,
    depends_on: Optional[List[Job]] = None,
    highmem_workers: bool = False,
    export_to_vcf: bool = False,
) -> List[Job]:
    """
    Add variant QC Hail-query jobs
    """
    vqsr_bucket = join(work_bucket, 'vqsr')
    
    depends_on = depends_on or []

    vqsred_vcf_path = join(vqsr_bucket, 'output.vcf.gz')
    if overwrite or not utils.file_exists(vqsred_vcf_path):
        vqsr_vcf_job = add_vqsr_jobs(
            b,
            combined_mt_path=raw_combined_mt_path,
            hard_filter_ht_path=hard_filter_ht_path,
            meta_ht_path=meta_ht_path,
            gvcf_count=sample_count,
            work_bucket=vqsr_bucket,
            web_bucket=join(web_bucket, 'vqsr'),
            vqsr_params_d=vqsr_params_d,
            scatter_count=scatter_count,
            output_vcf_path=vqsred_vcf_path,
            overwrite=overwrite,
            depends_on=depends_on,
        )
    else:
        vqsr_vcf_job = b.new_job('AS-VQSR [reuse]')

    job_name = 'AS-VQSR: load_vqsr'
    vqsr_split_ht_path = join(vqsr_bucket, 'vqsr-split.ht')
    if overwrite or not utils.file_exists(vqsr_split_ht_path):
        load_vqsr_job = add_job(
            b,
            f'{utils.SCRIPTS_DIR}/load_vqsr.py '
            f'--overwrite '
            f'--out-path {vqsr_split_ht_path} '
            f'--vqsr-vcf-path {vqsred_vcf_path} '
            f'--bucket {work_bucket} ',
            job_name=job_name,
            is_test=is_test,
            highmem=highmem_workers,
            num_workers=scatter_count,
            depends_on=depends_on,
        )
        load_vqsr_job.depends_on(vqsr_vcf_job)
    else:
        load_vqsr_job = b.new_job(f'{job_name} [reuse]')

    job_name = 'Var QC: generate QC annotations'
    allele_data_ht_path = join(work_bucket, 'allele-data.ht')
    # allele_data:
    #   nonsplit_alleles
    #   has_star
    #   variant_type
    #   n_alt_alleles
    #   allele_type
    #   was_mixed

    qc_ac_ht_path = join(work_bucket, 'qc-ac.ht')
    # ac_qc_samples_raw
    # ac_qc_samples_unrelated_raw
    # ac_release_samples_raw
    # ac_qc_samples_adj
    # ac_qc_samples_unrelated_adj
    # ac_release_samples_adj

    fam_stats_ht_path = join(work_bucket, 'fam-stats.ht') if ped_file else None

    if any(
        not utils.can_reuse(fp, overwrite)
        for fp in [allele_data_ht_path, qc_ac_ht_path]
        + ([fam_stats_ht_path] if fam_stats_ht_path else [])
    ):
        var_qc_anno_job = add_job(
            b,
            f'{utils.SCRIPTS_DIR}/generate_variant_qc_annotations.py '
            + f'{"--overwrite " if overwrite else ""}'
            + f'--mt {raw_combined_mt_path} '
            + f'--hard-filtered-samples-ht {hard_filter_ht_path} '
            + f'--meta-ht {meta_ht_path} '
            + f'--out-allele-data-ht {allele_data_ht_path} '
            + f'--out-qc-ac-ht {qc_ac_ht_path} '
            + (f'--out-fam-stats-ht {fam_stats_ht_path} ' if ped_file else '')
            + (f'--fam-file {ped_file} ' if ped_file else '')
            + f'--bucket {work_bucket} '
            + f'--n-partitions {scatter_count * 25}',
            job_name=job_name,
            is_test=is_test,
            highmem=highmem_workers,
            num_workers=scatter_count,
            depends_on=depends_on,
        )
        var_qc_anno_job.depends_on(*depends_on)
    else:
        var_qc_anno_job = b.new_job(f'{job_name} [reuse]')

    job_name = 'Var QC: generate frequencies'
    freq_ht_path = join(work_bucket, 'frequencies.ht')
    # InbreedingCoeff
    # freq
    # faf
    # popmax:
    #   AC
    #   AF
    #   AN
    #   homozygote_count
    #   pop
    #   faf95
    # qual_hists
    # raw_qual_hists

    if overwrite or not utils.file_exists(freq_ht_path):
        freq_job = add_job(
            b,
            f'{utils.SCRIPTS_DIR}/generate_freq_data.py --overwrite '
            f'--mt {raw_combined_mt_path} '
            f'--hard-filtered-samples-ht {hard_filter_ht_path} '
            f'--meta-ht {meta_ht_path} '
            f'--out-ht {freq_ht_path} '
            f'--bucket {work_bucket} ',
            job_name=job_name,
            is_test=is_test,
            highmem=highmem_workers,
            num_workers=scatter_count,
            depends_on=depends_on,
            long=True,
        )
        freq_job.depends_on(*depends_on)
    else:
        freq_job = b.new_job(f'{job_name} [reuse]')

    job_name = 'Making final MT'
    if not utils.can_reuse(out_filtered_combined_mt_path, overwrite):
        final_mt_j = add_job(
            b,
            f'{utils.SCRIPTS_DIR}/make_finalised_mt.py '
            f'--overwrite '
            f'--mt {raw_combined_mt_path} '
            f'--vqsr-ht {vqsr_split_ht_path} '
            f'--freq-ht {freq_ht_path} '
            f'--allele-data-ht {allele_data_ht_path} '
            f'--qc-ac-ht {qc_ac_ht_path} '
            f'--out-mt {out_filtered_combined_mt_path} '
            f'--meta-ht {meta_ht_path} ',
            job_name=job_name,
            is_test=is_test,
            num_workers=scatter_count,
            depends_on=depends_on,
        )
        final_mt_j.depends_on(load_vqsr_job, var_qc_anno_job, freq_job)
    else:
        final_mt_j = b.new_job(f'{job_name} [reuse]')

    jobs = [final_mt_j]
    if export_to_vcf:
        job_name = f'Making final VCF: prepare HT'
        logger.info(job_name)
        export_ht_path = join(work_bucket, 'export_vcf.ht')
        export_vcf_header_txt = join(work_bucket, 'export_vcf_header.txt')
        if not utils.can_reuse([export_ht_path, export_vcf_header_txt], overwrite):
            final_ht_j = add_job(
                b,
                f'{utils.SCRIPTS_DIR}/release_vcf_prepare_ht.py '
                f'--mt {out_filtered_combined_mt_path} '
                f'--out-ht {export_ht_path} '
                f'--out-vcf-header-txt {export_vcf_header_txt}',
                job_name=job_name,
                is_test=is_test,
                num_workers=scatter_count,
                depends_on=depends_on,
            )
        else:
            final_ht_j = b.new_job(f'{job_name} [reuse]')
        final_ht_j.depends_on(final_mt_j)
        jobs.append(final_mt_j)

        can_reuse_vcfs = False
        for chrom in list(map(str, range(1, 22 + 1))) + ['X', 'Y']:
            job_name = f'Making final VCF: HT to VCF for chr{chrom}'
            logger.info(job_name)
            vcf_path = out_filtered_vcf_ptrn_path.format(CHROM=chrom)
            if chrom == '1':
                can_reuse_vcfs = utils.can_reuse([vcf_path], overwrite)
            if not can_reuse_vcfs:
                j = add_job(
                    b,
                    f'{utils.SCRIPTS_DIR}/release_vcf_export_chrom.py '
                    f'--ht {export_ht_path} '
                    f'--vcf-header-txt {export_vcf_header_txt} '
                    f'--out-vcf {vcf_path} '
                    f'--name {project_name} '
                    f'--chromosome chr{chrom}',
                    job_name=job_name,
                    is_test=is_test,
                    num_workers=scatter_count,
                    depends_on=depends_on,
                )
            else:
                j = b.new_job(f'{job_name} [reuse]')
            jobs.append(j)
            j.depends_on(final_ht_j)
    return jobs
