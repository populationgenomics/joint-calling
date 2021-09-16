"""
Variant QC related hail-query jobs
"""

import uuid
from os.path import join
from typing import List, Optional, Dict, Tuple
import logging
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
    out_filtered_combined_mt_path: str,
    sample_count: int,
    ped_file: Optional[str],
    overwrite: bool,
    vqsr_params_d: Dict,
    scatter_count: int,
    is_test: bool,
    depends_on: Optional[List[Job]] = None,
    run_rf: bool = False,
) -> Job:
    """
    Add variant QC related hail-query jobs
    """
    max_age = '1h' if is_test else '12h'
    cluster1 = dataproc.setup_dataproc(
        b,
        max_age=max_age,
        packages=utils.DATAPROC_PACKAGES,
        num_secondary_workers=scatter_count,
        cluster_name='VarQC1',
        depends_on=depends_on,
    )

    rf_bucket = join(work_bucket, 'rf')
    vqsr_bucket = join(work_bucket, 'vqsr')

    fam_stats_ht_path = join(work_bucket, 'fam-stats.ht') if ped_file else None
    allele_data_ht_path = join(work_bucket, 'allele-data.ht')
    qc_ac_ht_path = join(work_bucket, 'qc-ac.ht')
    rf_result_ht_path = None

    job_name = 'Var QC: generate info'
    info_ht_path = join(work_bucket, 'info.ht')
    info_split_ht_path = join(work_bucket, 'info-split.ht')
    if any(
        not utils.can_reuse(fp, overwrite) for fp in [info_ht_path, info_split_ht_path]
    ):
        info_job = cluster1.add_job(
            f'{utils.SCRIPTS_DIR}/generate_info_ht.py --overwrite '
            f'--mt {raw_combined_mt_path} '
            f'--out-info-ht {info_ht_path} '
            f'--out-split-info-ht {info_split_ht_path}',
            job_name=job_name,
        )
        if depends_on:
            info_job.depends_on(*depends_on)
    else:
        info_job = b.new_job(f'{job_name} [reuse]')

    max_age = '1h' if is_test else '12h'
    cluster2 = dataproc.setup_dataproc(
        b,
        max_age=max_age,
        packages=utils.DATAPROC_PACKAGES,
        num_secondary_workers=scatter_count,
        cluster_name='VarQC2',
        depends_on=depends_on,
    )

    job_name = 'Var QC: generate annotations'
    if any(
        not utils.can_reuse(fp, overwrite)
        for fp in [allele_data_ht_path, qc_ac_ht_path]
    ):
        var_qc_anno_job = cluster2.add_job(
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
        )
        if depends_on:
            var_qc_anno_job.depends_on(*depends_on)
    else:
        var_qc_anno_job = b.new_job(f'{job_name} [reuse]')

    max_age = '2h' if is_test else '24h'
    cluster3 = dataproc.setup_dataproc(
        b,
        max_age=max_age,
        packages=utils.DATAPROC_PACKAGES,
        num_secondary_workers=scatter_count,
        cluster_name='VarQC3',
        depends_on=depends_on,
    )

    job_name = 'Var QC: generate frequencies'
    freq_ht_path = join(work_bucket, 'frequencies.ht')
    if overwrite or not utils.file_exists(freq_ht_path):
        freq_job = cluster3.add_job(
            f'{utils.SCRIPTS_DIR}/generate_freq_data.py --overwrite '
            f'--mt {raw_combined_mt_path} '
            f'--hard-filtered-samples-ht {hard_filter_ht_path} '
            f'--meta-ht {meta_ht_path} '
            f'--out-ht {freq_ht_path} '
            f'--bucket {work_bucket} ',
            job_name=job_name,
        )
        if depends_on:
            freq_job.depends_on(*depends_on)
    else:
        freq_job = b.new_job('{job_name} [reuse]')

    job_name = 'Var QC: create RF annotations'
    rf_annotations_ht_path = join(work_bucket, 'rf-annotations.ht')
    if overwrite or not utils.file_exists(rf_annotations_ht_path):
        rf_anno_job = cluster3.add_job(
            f'{utils.SCRIPTS_DIR}/create_rf_annotations.py --overwrite '
            f'--info-split-ht {info_split_ht_path} '
            f'--freq-ht {freq_ht_path} '
            + (f'--fam-stats-ht {fam_stats_ht_path} ' if fam_stats_ht_path else '')
            + f'--allele-data-ht {allele_data_ht_path} '
            f'--qc-ac-ht {qc_ac_ht_path} '
            f'--bucket {work_bucket} '
            f'--use-adj-genotypes '
            f'--out-ht {rf_annotations_ht_path} '
            + f'--n-partitions {scatter_count * 25}',
            job_name=job_name,
        )
        rf_anno_job.depends_on(freq_job, var_qc_anno_job, info_job)
    else:
        rf_anno_job = b.new_job(f'{job_name} [reuse]')

    if run_rf:
        job_name = 'Random forest'
        rf_result_ht_path = join(work_bucket, 'rf-result.ht')
        rf_model_id = f'rf_{str(uuid.uuid4())[:8]}'
        if overwrite or not utils.file_exists(rf_result_ht_path):
            rf_job = cluster3.add_job(
                f'{utils.SCRIPTS_DIR}/random_forest.py --overwrite '
                f'--annotations-ht {rf_annotations_ht_path} '
                f'--bucket {work_bucket} '
                f'--use-adj-genotypes '
                f'--out-results-ht {rf_result_ht_path} '
                f'--out-model-id {rf_model_id} ',
                job_name=job_name,
            )
            rf_job.depends_on(rf_anno_job)
        else:
            rf_job = b.new_job(f'{job_name} [reuse]')

        eval_job, final_filter_ht_path = make_rf_eval_jobs(
            b=b,
            cluster=cluster3,
            combined_mt_path=raw_combined_mt_path,
            info_split_ht_path=info_split_ht_path,
            rf_result_ht_path=rf_result_ht_path,
            rf_annotations_ht_path=rf_annotations_ht_path,
            fam_stats_ht_path=fam_stats_ht_path,
            freq_ht_path=freq_ht_path,
            rf_model_id=rf_model_id,
            work_bucket=rf_bucket,
            overwrite=overwrite,
            depends_on=[rf_job, freq_job, rf_anno_job],
        )

    else:
        vqsred_vcf_path = join(vqsr_bucket, 'output.vcf.gz')
        if overwrite or not utils.file_exists(vqsred_vcf_path):
            vqsr_vcf_job = make_vqsr_jobs(
                b,
                combined_mt_path=raw_combined_mt_path,
                hard_filter_ht_path=hard_filter_ht_path,
                meta_ht_path=meta_ht_path,
                gvcf_count=sample_count,
                work_bucket=vqsr_bucket,
                web_bucket=join(web_bucket, 'vqsr'),
                depends_on=depends_on or [],
                vqsr_params_d=vqsr_params_d,
                scatter_count=scatter_count,
                output_vcf_path=vqsred_vcf_path,
                overwrite=overwrite,
            )
        else:
            vqsr_vcf_job = b.new_job('AS-VQSR [reuse]')

        max_age = '1h' if is_test else '12h'
        cluster3 = dataproc.setup_dataproc(
            b,
            max_age=max_age,
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=scatter_count,
            cluster_name='VQeval',
            depends_on=[vqsr_vcf_job],
        )
        final_filter_ht_path = join(vqsr_bucket, 'final-filter.ht')
        eval_job = make_vqsr_eval_jobs(
            b=b,
            cluster=cluster3,
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
            final_gathered_vcf_job=vqsr_vcf_job,
            rf_anno_job=rf_anno_job,
            output_ht_path=final_filter_ht_path,
        )
        eval_job.depends_on(vqsr_vcf_job, rf_anno_job, info_job)

    job_name = 'Making final MT'
    if not utils.can_reuse(out_filtered_combined_mt_path, overwrite):
        final_job = cluster3.add_job(
            f'{utils.SCRIPTS_DIR}/make_finalised_mt.py --overwrite '
            f'--mt {raw_combined_mt_path} '
            f'--final-filter-ht {final_filter_ht_path} '
            f'--out-mt {out_filtered_combined_mt_path} '
            f'--meta-ht {meta_ht_path} ',
            job_name=job_name,
        )
        final_job.depends_on(eval_job)
    else:
        final_job = b.new_job(f'{job_name} [reuse]')

    return final_job


def make_rf_eval_jobs(
    b: hb.Batch,
    cluster,
    combined_mt_path: str,
    info_split_ht_path: str,
    rf_result_ht_path: str,
    rf_annotations_ht_path: str,
    fam_stats_ht_path: Optional[str],
    freq_ht_path: str,
    rf_model_id: str,
    work_bucket: str,
    overwrite: bool,
    depends_on: Optional[List[Job]],
) -> Tuple[Job, str]:
    """
    Make jobs that do evaluation RF model and applies the final filters

    Returns the final_filter Job object and the path to the final filter HT
    """
    job_name = 'RF: evaluation'
    score_bin_ht_path = join(work_bucket, 'rf-score-bin.ht')
    score_bin_agg_ht_path = join(work_bucket, 'rf-score-agg-bin.ht')
    if overwrite or not utils.file_exists(score_bin_ht_path):
        eval_job = cluster.add_job(
            b,
            f'{utils.SCRIPTS_DIR}/evaluation.py --overwrite '
            f'--mt {combined_mt_path} '
            f'--rf-annotations-ht {rf_annotations_ht_path} '
            f'--info-split-ht {info_split_ht_path} '
            + (f'--fam-stats-ht {fam_stats_ht_path} ' if fam_stats_ht_path else '')
            + f'--rf-results-ht {rf_result_ht_path} '
            f'--bucket {work_bucket} '
            f'--out-bin-ht {score_bin_ht_path} '
            f'--out-aggregated-bin-ht {score_bin_agg_ht_path} '
            f'--run-sanity-checks ',
            job_name=job_name,
        )
        if depends_on:
            eval_job.depends_on(*depends_on)
    else:
        eval_job = b.new_job(f'{job_name} [reuse]')

    job_name = 'RF: final filter'
    final_filter_ht_path = join(work_bucket, 'final-filter.ht')
    if overwrite or not utils.file_exists(final_filter_ht_path):
        final_filter_job = cluster.add_job(
            f'{utils.SCRIPTS_DIR}/final_filter.py --overwrite '
            f'--out-final-filter-ht {final_filter_ht_path} '
            f'--model-id {rf_model_id} '
            f'--model-name RF '
            f'--score-name RF '
            f'--info-split-ht {info_split_ht_path} '
            f'--freq-ht {freq_ht_path} '
            f'--score-bin-ht {score_bin_ht_path} '
            f'--score-bin-agg-ht {score_bin_agg_ht_path} ' + f'--bucket {work_bucket} ',
            job_name=job_name,
        )
        final_filter_job.depends_on(eval_job)
    else:
        final_filter_job = b.new_job(f'{job_name} [reuse]')

    return final_filter_job, final_filter_ht_path


def make_vqsr_eval_jobs(
    b: hb.Batch,
    cluster,
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
    final_gathered_vcf_job: Job,
    rf_anno_job: Job,
    output_ht_path: str,
) -> Job:
    """
    Make jobs that do evaluation VQSR model and applies the final filters

    Returns the final_filter Job object and the path to the final filter HT
    """
    job_name = 'AS-VQSR: load_vqsr'
    vqsr_filters_split_ht_path = join(work_bucket, 'vqsr-filters-split.ht')
    if overwrite or not utils.file_exists(vqsr_filters_split_ht_path):
        load_vqsr_job = cluster.add_job(
            f'{utils.SCRIPTS_DIR}/load_vqsr.py --overwrite '
            f'--split-multiallelic '
            f'--out-path {vqsr_filters_split_ht_path} '
            f'--vqsr-vcf-path {final_gathered_vcf_path} '
            f'--bucket {work_bucket} ',
            job_name=job_name,
        )
        load_vqsr_job.depends_on(final_gathered_vcf_job)
    else:
        load_vqsr_job = b.new_job(f'{job_name} [reuse]')

    job_name = 'AS-VQSR: evaluation'
    score_bin_ht_path = join(work_bucket, 'vqsr-score-bin.ht')
    score_bin_agg_ht_path = join(work_bucket, 'vqsr-score-agg-bin.ht')
    if (
        overwrite
        or not utils.file_exists(score_bin_ht_path)
        or not utils.file_exists(score_bin_agg_ht_path)
    ):
        eval_job = cluster.add_job(
            f'{utils.SCRIPTS_DIR}/evaluation.py --overwrite '
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
            job_name=job_name,
        )
        eval_job.depends_on(load_vqsr_job, rf_anno_job)
    else:
        eval_job = b.new_job(f'{job_name} [reuse]')

    job_name = 'AS-VQSR: final filter'
    vqsr_model_id = 'vqsr_model'
    if not utils.file_exists(output_ht_path):
        final_filter_job = cluster.add_job(
            f'{utils.SCRIPTS_DIR}/final_filter.py --overwrite '
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
            job_name=job_name,
        )
        final_filter_job.depends_on(eval_job)
    else:
        final_filter_job = b.new_job(f'{job_name} [reuse]')
    return final_filter_job
