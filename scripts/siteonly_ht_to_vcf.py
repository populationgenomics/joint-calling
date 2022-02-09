#!/usr/bin/env python

"""
Convert annotated site-only table to a sites-only VCF.
Essentially a verbatim copy of: hail-ukbb-200k-callset:mt_to_vcf.py
"""

import logging
import click
import hail as hl
from gnomad.utils.vcf import adjust_vcf_incompatible_types
from gnomad.utils.sparse_mt import default_compute_info
from joint_calling import _version
from joint_calling.utils import get_validation_callback, init_hail, file_exists
from joint_calling import utils

logger = logging.getLogger('vqsr_qc')
logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--ht',
    'ht_path',
    required=True,
    callback=get_validation_callback(ext='ht', must_exist=True),
    help='path to annotated Hail table',
)
@click.option(
    '-o',
    'output_path',
    required=True,
)
@click.option(
    '--local-tmp-dir',
    'local_tmp_dir',
    help='local directory for temporary files and Hail logs (must be local).',
)
@click.option(
    '--overwrite/--reuse',
    'overwrite',
    is_flag=True,
    help='if an intermediate or a final file exists, skip running the code '
    'that generates it.',
)
@click.option(
    '--hail-billing',
    'hail_billing',
    help='Hail billing account ID.',
)
def main(
    ht_path: str,
    output_path: str,
    local_tmp_dir: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):  # pylint: disable=missing-function-docstring
    init_hail(__file__, local_tmp_dir)
    
    if utils.can_reuse(output_path, overwrite):
        return
    
    logger.info(
        f'Exporting sites-only VCF to {output_path} to run in the VQSR pipeline'
    )
    hl.export_vcf(hl.read_table(ht_path), output_path)
    logger.info('Successfully exported sites-only VCF')


if __name__ == '__main__':
    main()  # pylint: disable=E1120
