"""
Run the evaluation script
"""

import os
from os.path import join, dirname, abspath
import logging
import hailtop.batch as hb
from analysis_runner import dataproc

from joint_calling import utils

logger = logging.getLogger('test-joint-calling')
logger.setLevel('INFO')


def main():
    """
    Drive a Hail Batch workflow that creates and submits jobs. A job usually runs
    either a Hail Query script from the scripts folder in this repo using a Dataproc
    cluster; or a GATK command using the GATK or Gnarly image.
    """
    backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.environ.get('HAIL_BUCKET').replace('gs://', ''),
    )
    b = hb.Batch('Joint Calling', backend=backend)
    scripts_dir = abspath(join(dirname(__file__), '..', 'scripts'))

    raw_combined_mt_path = 'gs://cpg-tob-wgs-temporary/joint_vcf/v1/raw/genomes.mt/'
    info_split_ht_path = (
        'gs://cpg-tob-wgs-temporary/joint_vcf/v1/work/variant_qc/info.ht/'
    )
    freq_ht_path = (
        'gs://cpg-tob-wgs-temporary/joint_vcf/v1/work/variant_qc/frequencies.ht/'
    )
    rf_model_id = 'rf_375c9463'
    rf_result_ht_path = (
        'gs://cpg-tob-wgs-temporary/joint_vcf/v1/work/variant_qc/rf_result.ht/'
    )
    rf_annotations_ht_path = (
        'gs://cpg-tob-wgs-temporary/joint_vcf/v1/work/variant_qc/rf_annotations-adj.ht/'
    )

    eval_bucket = (
        'gs://cpg-tob-wgs-temporary/joint_vcf/v1/work/variant_qc/evaluation/rf/'
    )
    score_bin_ht_path = join(eval_bucket, 'score-bin.ht')
    score_bin_agg_ht_path = join(eval_bucket, 'score-agg-bin.ht')
    fam_stats_ht_path = (
        'gs://cpg-tob-wgs-temporary/joint_vcf/v1/work/variant_qc/fam-stats.ht'
    )

    if not utils.file_exists(score_bin_ht_path) or not utils.file_exists(
        score_bin_agg_ht_path
    ):
        eval_job = dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/evaluation.py --reuse '
            f'--info-split-ht {info_split_ht_path} '
            f'--rf-results-ht {rf_result_ht_path} '
            f'--rf-annotations-ht {rf_annotations_ht_path} '
            f'--fam-stats-ht {fam_stats_ht_path} '
            f'--mt {raw_combined_mt_path} '
            f'--bucket {eval_bucket} '
            f'--out-bin-ht {score_bin_ht_path} '
            f'--out-aggregated-bin-ht {score_bin_agg_ht_path} '
            f'--run-sanity-checks ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            depends_on=[],
            job_name='RF: evaluation',
        )
    else:
        eval_job = b.new_job('RF: evaluation')

    final_filter_ht_path = join(eval_bucket, 'final-filter.ht')
    if not utils.file_exists(final_filter_ht_path):
        dataproc.hail_dataproc_job(
            b,
            f'{scripts_dir}/final_filter.py --overwrite '
            f'--out-final-filter-ht {final_filter_ht_path} '
            f'--model-id {rf_model_id} '
            f'--model-name RF '
            f'--score-name RF '
            f'--info-split-ht {info_split_ht_path} '
            f'--freq-ht {freq_ht_path} '
            f'--score-bin-ht {score_bin_ht_path} '
            f'--score-bin-agg-ht {score_bin_agg_ht_path} '
            f'--bucket {eval_bucket} ',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            depends_on=[eval_job],
            job_name='RF: final filter',
        )

    b.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
