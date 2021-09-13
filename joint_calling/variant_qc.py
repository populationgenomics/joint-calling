"""
Variant QC related hail-query jobs
"""

import uuid
from os.path import join
from typing import List, Optional, Dict, Tuple
import logging
import pandas as pd
import hailtop.batch as hb
from hailtop.batch.job import Job

from analysis_runner import dataproc
from joint_calling import utils
from joint_calling.vqsr import make_vqsr_jobs

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
    samples_df: pd.DataFrame,
    scripts_dir: str,
    ped_file: Optional[str],
    overwrite: bool,
    vqsr_params_d: Dict,
    scatter_count: int,
    depends_on: Optional[List[Job]] = None,
    run_rf: bool = False,
) -> Tuple[Job, str]:
    """
    Add variant QC related hail-query jobs
    """
    rf_bucket = join(work_bucket, 'rf')
    vqsr_bucket = join(work_bucket, 'vqsr')

    fam_stats_ht_path = join(work_bucket, 'fam-stats.ht') if ped_file else None
    allele_data_ht_path = join(work_bucket, 'allele-data.ht')
    qc_ac_ht_path = join(work_bucket, 'qc-ac.ht')
    rf_result_ht_path = None

    info_ht_path = join(work_bucket, 'info.ht')
    info_split_ht_path = join(work_bucket, 'info-split.ht')
    if any(
        not utils.can_reuse(fp, overwrite) for fp in [info_ht_path, info_split_ht_path]
    ):
        info_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/generate_info_ht.py --overwrite '
            f'--mt {raw_combined_mt_path} '
            f'--out-info-ht {info_ht_path} '
            f'--out-split-info-ht {info_split_ht_path}',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=depends_on or [],
            job_name='Var QC: generate info',
        )
    else:
        info_job = b.new_job('Var QC: generate info [reuse]')

    if any(
        not utils.can_reuse(fp, overwrite)
        for fp in [allele_data_ht_path, qc_ac_ht_path]
    ):
        var_qc_anno_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/generate_variant_qc_annotations.py '
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
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=depends_on or [],
            job_name='Var QC: generate annotations',
            vep='GRCh38',
        )
    else:
        var_qc_anno_job = b.new_job('Var QC: generate annotations [reuse]')

    freq_ht_path = join(work_bucket, 'frequencies.ht')
    if overwrite or not utils.file_exists(freq_ht_path):
        freq_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/generate_freq_data.py --overwrite '
            f'--mt {raw_combined_mt_path} '
            f'--hard-filtered-samples-ht {hard_filter_ht_path} '
            f'--meta-ht {meta_ht_path} '
            f'--out-ht {freq_ht_path} '
            f'--bucket {work_bucket} ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=depends_on or [],
            job_name='Var QC: generate frequencies',
        )
    else:
        freq_job = b.new_job('Var QC: generate frequencies [reuse]')

    rf_annotations_ht_path = join(work_bucket, 'rf-annotations.ht')
    if overwrite or not utils.file_exists(rf_annotations_ht_path):
        rf_anno_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/create_rf_annotations.py --overwrite '
            f'--info-split-ht {info_split_ht_path} '
            f'--freq-ht {freq_ht_path} '
            + (f'--fam-stats-ht {fam_stats_ht_path} ' if fam_stats_ht_path else '')
            + f'--allele-data-ht {allele_data_ht_path} '
            f'--qc-ac-ht {qc_ac_ht_path} '
            f'--bucket {work_bucket} '
            f'--use-adj-genotypes '
            f'--out-ht {rf_annotations_ht_path} '
            + f'--n-partitions {scatter_count * 25}',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=[freq_job, var_qc_anno_job, info_job],
            job_name='Var QC: create RF annotations',
        )
    else:
        rf_anno_job = b.new_job('Var QC: create RF annotations [reuse]')

    if run_rf:
        rf_result_ht_path = join(work_bucket, 'rf-result.ht')
        rf_model_id = f'rf_{str(uuid.uuid4())[:8]}'
        if overwrite or not utils.file_exists(rf_result_ht_path):
            rf_job = dataproc.hail_dataproc_job(
                b,
                f'{scripts_dir}/random_forest.py --overwrite '
                f'--annotations-ht {rf_annotations_ht_path} '
                f'--bucket {work_bucket} '
                f'--use-adj-genotypes '
                f'--out-results-ht {rf_result_ht_path} '
                f'--out-model-id {rf_model_id} ',
                max_age='8h',
                packages=utils.DATAPROC_PACKAGES,
                num_secondary_workers=scatter_count,
                depends_on=[rf_anno_job],
                job_name='Random forest',
            )
        else:
            rf_job = b.new_job('Random forest [reuse]')

        final_job, final_filter_ht_path = make_rf_eval_jobs(
            b=b,
            combined_mt_path=raw_combined_mt_path,
            info_split_ht_path=info_split_ht_path,
            rf_result_ht_path=rf_result_ht_path,
            rf_annotations_ht_path=rf_annotations_ht_path,
            fam_stats_ht_path=fam_stats_ht_path,
            freq_ht_path=freq_ht_path,
            rf_model_id=rf_model_id,
            work_bucket=rf_bucket,
            overwrite=overwrite,
            scripts_dir=scripts_dir,
            depends_on=[rf_job],
            scatter_count=scatter_count,
        )

    else:
        vqsred_vcf_path = join(vqsr_bucket, 'output.vcf.gz')
        if overwrite or not utils.file_exists(vqsred_vcf_path):
            final_gathered_vcf_job = make_vqsr_jobs(
                b,
                combined_mt_path=raw_combined_mt_path,
                hard_filter_ht_path=hard_filter_ht_path,
                meta_ht_path=meta_ht_path,
                gvcf_count=len(samples_df),
                work_bucket=vqsr_bucket,
                web_bucket=join(web_bucket, 'vqsr'),
                depends_on=depends_on or [],
                scripts_dir=scripts_dir,
                vqsr_params_d=vqsr_params_d,
                scatter_count=scatter_count,
                output_vcf_path=vqsred_vcf_path,
                overwrite=overwrite,
            )
        else:
            final_gathered_vcf_job = b.new_job('VQSR [reuse]')

        final_filter_ht_path = join(vqsr_bucket, 'final-filter.ht')
        final_job = make_vqsr_eval_jobs(
            b=b,
            combined_mt_path=raw_combined_mt_path,
            rf_annotations_ht_path=rf_annotations_ht_path,
            info_split_ht_path=info_split_ht_path,
            final_gathered_vcf_path=vqsred_vcf_path,
            rf_result_ht_path=rf_result_ht_path,
            fam_stats_ht_path=fam_stats_ht_path,
            freq_ht_path=freq_ht_path,
            work_bucket=vqsr_bucket,
            analysis_bucket=join(web_bucket, 'vqsr'),
            overwrite=overwrite,
            scripts_dir=scripts_dir,
            final_gathered_vcf_job=final_gathered_vcf_job,
            rf_anno_job=rf_anno_job,
            scatter_count=scatter_count,
            output_ht_path=final_filter_ht_path,
        )
    return final_job, final_filter_ht_path


def make_rf_eval_jobs(
    b: hb.Batch,
    combined_mt_path: str,
    info_split_ht_path: str,
    rf_result_ht_path: str,
    rf_annotations_ht_path: str,
    fam_stats_ht_path: Optional[str],
    freq_ht_path: str,
    rf_model_id: str,
    work_bucket: str,
    overwrite: bool,
    scripts_dir: str,
    depends_on: Optional[List[Job]],
    scatter_count: int,
) -> Tuple[Job, str]:
    """
    Make jobs that do evaluation RF model and applies the final filters

    Returns the final_filter Job object and the path to the final filter HT
    """
    score_bin_ht_path = join(work_bucket, 'rf-score-bin.ht')
    score_bin_agg_ht_path = join(work_bucket, 'rf-score-agg-bin.ht')
    if overwrite or not utils.file_exists(score_bin_ht_path):
        eval_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/evaluation.py --overwrite '
            f'--mt {combined_mt_path} '
            f'--rf-annotations-ht {rf_annotations_ht_path} '
            f'--info-split-ht {info_split_ht_path} '
            + (f'--fam-stats-ht {fam_stats_ht_path} ' if fam_stats_ht_path else '')
            + f'--rf-results-ht {rf_result_ht_path} '
            f'--bucket {work_bucket} '
            f'--out-bin-ht {score_bin_ht_path} '
            f'--out-aggregated-bin-ht {score_bin_agg_ht_path} '
            f'--run-sanity-checks ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=depends_on,
            job_name='RF: evaluation',
        )
    else:
        eval_job = b.new_job('RF: evaluation [reuse]')

    final_filter_ht_path = join(work_bucket, 'final-filter.ht')
    if overwrite or not utils.file_exists(final_filter_ht_path):
        final_filter_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/final_filter.py --overwrite '
            f'--out-final-filter-ht {final_filter_ht_path} '
            f'--model-id {rf_model_id} '
            f'--model-name RF '
            f'--score-name RF '
            f'--info-split-ht {info_split_ht_path} '
            f'--freq-ht {freq_ht_path} '
            f'--score-bin-ht {score_bin_ht_path} '
            f'--score-bin-agg-ht {score_bin_agg_ht_path} ' + f'--bucket {work_bucket} ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=[eval_job],
            job_name='RF: final filter',
        )
    else:
        final_filter_job = b.new_job('RF: final filter [reuse]')

    return final_filter_job, final_filter_ht_path


def make_vqsr_eval_jobs(
    b: hb.Batch,
    combined_mt_path: str,
    rf_annotations_ht_path: str,
    info_split_ht_path: str,
    final_gathered_vcf_path: str,
    rf_result_ht_path: Optional[str],
    fam_stats_ht_path: Optional[str],
    freq_ht_path: str,
    work_bucket: str,
    analysis_bucket: str,  # pylint: disable=unused-argument
    overwrite: bool,
    scripts_dir: str,
    final_gathered_vcf_job: Job,
    rf_anno_job: Job,
    scatter_count: int,
    output_ht_path: str,
) -> Tuple[Job, str]:
    """
    Make jobs that do evaluation VQSR model and applies the final filters

    Returns the final_filter Job object and the path to the final filter HT
    """
    vqsr_filters_split_ht_path = join(work_bucket, 'vqsr-filters-split.ht')
    if overwrite or not utils.file_exists(vqsr_filters_split_ht_path):
        load_vqsr_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/load_vqsr.py --overwrite '
            f'--split-multiallelic '
            f'--out-path {vqsr_filters_split_ht_path} '
            f'--vqsr-vcf-path {final_gathered_vcf_path} '
            f'--bucket {work_bucket} ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=[final_gathered_vcf_job],
            job_name='VQSR: load_vqsr',
        )
    else:
        load_vqsr_job = b.new_job('VQSR: load_vqsr [reuse]')

    score_bin_ht_path = join(work_bucket, 'vqsr-score-bin.ht')
    score_bin_agg_ht_path = join(work_bucket, 'vqsr-score-agg-bin.ht')
    if (
        overwrite
        or not utils.file_exists(score_bin_ht_path)
        or not utils.file_exists(score_bin_agg_ht_path)
    ):
        eval_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/evaluation.py --overwrite '
            f'--mt {combined_mt_path} '
            f'--rf-annotations-ht {rf_annotations_ht_path} '
            f'--info-split-ht {info_split_ht_path} '
            + (f'--fam-stats-ht {fam_stats_ht_path} ' if fam_stats_ht_path else '')
            + (
                f'--rf-result-ht {rf_result_ht_path} '
                if (rf_annotations_ht_path and rf_result_ht_path)
                else ''
            )
            + f'--vqsr-filters-split-ht {vqsr_filters_split_ht_path} '
            f'--bucket {work_bucket} '
            f'--out-bin-ht {score_bin_ht_path} '
            f'--out-aggregated-bin-ht {score_bin_agg_ht_path} '
            f'--run-sanity-checks ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=[load_vqsr_job, rf_anno_job],
            job_name='VQSR: evaluation',
        )
    else:
        eval_job = b.new_job('VQSR: evaluation [reuse]')

    vqsr_model_id = 'vqsr_model'
    if not utils.file_exists(output_ht_path):
        final_filter_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/final_filter.py --overwrite '
            f'--out-final-filter-ht {output_ht_path} '
            f'--vqsr-filters-split-ht {vqsr_filters_split_ht_path} '
            f'--model-id {vqsr_model_id} '
            f'--model-name VQSR '
            f'--score-name AS_VQSLOD '
            f'--info-split-ht {info_split_ht_path} '
            f'--freq-ht {freq_ht_path} '
            f'--score-bin-ht {score_bin_ht_path} '
            f'--score-bin-agg-ht {score_bin_agg_ht_path} '
            f'--bucket {work_bucket} ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            depends_on=[eval_job],
            job_name='VQSR: final filter',
        )
    else:
        final_filter_job = b.new_job('VQSR: final filter [reuse]')
    return final_filter_job
