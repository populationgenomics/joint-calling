"""
Hail Batch workflow to perform joint calling, sample QC, and variant QC with VQSR and 
random forest methods on a WGS germline callset.

The workflow is parametrised by the access level, the dataset name, 
batch names and the output version.

It must be only run with the CPG analysis-runner:
https://github.com/populationgenomics/analysis-runner (see helper script `driver_for_analysis_runner.sh` for analysis-runner submissions)
"""

import os
from os.path import join, dirname, abspath
import logging
import hailtop.batch as hb
from analysis_runner import dataproc
from joint_calling import utils

logger = logging.getLogger('joint-calling')
logger.setLevel('INFO')


def main():  # pylint: disable=too-many-arguments,too-many-locals,too-many-statements
    """
    Drive a Hail Batch workflow that creates and submits jobs. A job usually runs
    either a Hail Query script from the scripts folder in this repo using a Dataproc
    cluster; or a GATK command using the GATK or Gnarly image.
    """
    logger.info(f'Enable random forest')
    logger.info(f'Enable VQSR')

    billing_project = os.getenv('HAIL_BILLING_PROJECT')
    hail_bucket = os.environ.get('HAIL_BUCKET')
    logger.info(
        f'Starting hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=hail_bucket.replace('gs://', ''),
    )
    b = hb.Batch(
        f'Test mt2vcf',
        backend=backend,
    )
    scripts_dir = abspath(join(dirname(__file__), '..', 'scripts'))

    mt_path = 'gs://cpg-tob-wgs-test/mt/v2-raw.mt'
    meta_ht_path = 'gs://cpg-tob-wgs-main-metadata/joint-calling/v2/meta.ht/'
    hard_filter_ht_path = (
        'gs://cpg-tob-wgs-test-tmp/joint-calling/v3/sample_qc/hard_filters.ht/'
    )
    combined_vcf_path = (
        'gs://cpg-tob-wgs-test-tmp/joint-calling/v2/variant_qc/vqsr/input-test.vcf.gz'
    )
    dataproc.hail_dataproc_job(
        b,
        f'{scripts_dir}/mt_to_vcf.py --overwrite '
        f'--mt {mt_path} '
        f'--meta-ht {meta_ht_path} '
        f'--hard-filtered-samples-ht {hard_filter_ht_path} '
        f'-o {combined_vcf_path} ',
        max_age='8h',
        packages=utils.DATAPROC_PACKAGES,
        num_secondary_workers=50,
        secondary_worker_boot_disk_size=10,
        job_name='MT to VCF',
    )

    b.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
