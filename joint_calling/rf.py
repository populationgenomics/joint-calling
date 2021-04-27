"""
Create jobs to create and apply a Random Forest model
"""

from os.path import join
from typing import List, Optional
import logging
import hailtop.batch as hb
from hailtop.batch.job import Job
from analysis_runner import dataproc

from joint_calling import utils

logger = logging.getLogger('joint-calling')
logger.setLevel('INFO')


def make_rf_jobs(
    b: hb.Batch,
    combined_mt_path: str,
    hard_filtered_samples_ht_path: str,
    meta_ht_path: str,
    work_bucket: str,
    depends_on: Optional[List[Job]],
) -> Job:
    """
    :param b: Batch object to add jobs to
    :param combined_mt_path: path to a Matrix Table combined with the Hail VCF combiner
    :param hard_filtered_samples_ht_path: a Table containing only samples that failed the hard filters
    :param meta_ht_path: a Table with meta data from sample_qc.py
    :param work_bucket: bucket for intermediate files
    :param depends_on: job that the created jobs should only run after
    :return: a Job object of the last created job
    """
    freq_ht_path = join(work_bucket, 'frequencies.ht')
    info_ht_path = join(work_bucket, 'info.ht')
    allele_data_ht_path = join(work_bucket, 'allele_data.ht')
    qc_ac_ht_path = join(work_bucket, 'qc_ac.ht')
    rf_result_ht_path = join(work_bucket, 'rf_result.ht')
    if not all(
        utils.file_exists(fp)
        for fp in [info_ht_path, allele_data_ht_path, qc_ac_ht_path]
    ):
        rf_anno_job = dataproc.hail_dataproc_job(
            b,
            f'scripts/run_python_script.py '
            f'generate_qc_annotations.py --overwrite '
            f'--split-multiallelic '
            f'--mt {combined_mt_path} '
            f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
            f'--meta-ht {meta_ht_path} '
            f'--out-info-ht {info_ht_path} '
            f'--out-allele-data-ht {allele_data_ht_path} '
            f'--out-qc-ac-ht {qc_ac_ht_path} '
            f'--bucket {work_bucket} ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            depends_on=depends_on,
            job_name='RF: gen QC anno',
            vep='GRCh38',
        )
    else:
        rf_anno_job = b.new_job('RF: gen QC anno')

    if not utils.file_exists(freq_ht_path):
        rf_freq_data_job = dataproc.hail_dataproc_job(
            b,
            f'scripts/run_python_script.py '
            f'generate_freq_data.py --overwrite '
            f'--mt {combined_mt_path} '
            f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
            f'--meta-ht {meta_ht_path} '
            f'--out-ht {freq_ht_path} '
            f'--bucket {work_bucket} ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            depends_on=depends_on,
            job_name='RF: gen freq data',
        )
    else:
        rf_freq_data_job = b.new_job('RF: gen freq data')

    if not utils.file_exists(rf_result_ht_path):
        rf_job = dataproc.hail_dataproc_job(
            b,
            f'scripts/run_python_script.py '
            f'random_forest.py --overwrite '
            f'--info-ht {info_ht_path} '
            f'--freq-ht {freq_ht_path} '
            f'--allele-data-ht {allele_data_ht_path} '
            f'--qc-ac-ht {qc_ac_ht_path} '
            f'--bucket {work_bucket} '
            '--use-adj-genotypes '
            f'--out-ht {rf_result_ht_path} ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            depends_on=[rf_freq_data_job, rf_anno_job],
            job_name='RF: run',
        )
    else:
        rf_job = b.new_job('RF: run')

    return rf_job
