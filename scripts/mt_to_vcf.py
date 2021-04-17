#!/usr/bin/env python

"""
This script takes a path to a matrix table, and then:
 - exports a sites-only VCF using Hail Query,
 - run a pipeline on Cromwell
 - Collect the results
"""

import logging
import click
import hail as hl

from joint_calling import _version
from joint_calling.utils import get_validation_callback, init_hail, file_exists
from joint_calling.mt_to_vcf import mt_to_sites_only_mt

logger = logging.getLogger('vqsr_qc')
logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--mt',
    'mt_path',
    required=True,
    callback=get_validation_callback(ext='mt', must_exist=True),
    help='path to the input MatrixTable',
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
    '--overwrite',
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
@click.option(
    '--partitions',
    'partitions',
    default=5000,
    help='Number of partitions for Hail distributed computing',
)
def main(
    mt_path: str,
    output_path: str,
    local_tmp_dir: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
    partitions: int,
):
    """
    Expects hail service to already be initialised
    """
    init_hail('variant_qc', local_tmp_dir)

    logger.info(f'Loading matrix table from "{mt_path}"')
    mt = hl.read_matrix_table(mt_path).key_rows_by('locus', 'alleles')

    if file_exists(output_path):
        if overwrite:
            logger.info(f'Output file {output_path} exists and will be overwritten')
        else:
            logger.info(
                f'Output file {output_path} exists, use --overwrite to overwrite'
            )
            return
    export_sites_only_vcf(mt=mt, output_path=output_path, partitions=partitions)


def export_sites_only_vcf(mt: hl.MatrixTable, output_path: str, partitions: int = 5000):
    """
    Take initial matrix table, convert to sites-only matrix table, then export to vcf
    """
    logger.info('Converting matrix table to sites-only matrix table')
    final_mt = mt_to_sites_only_mt(mt, partitions)

    # export vcf, and return the path

    logger.info(
        f"Exporting sites-only VCF to '{output_path}' to run in the VQSR pipeline"
    )
    hl.export_vcf(final_mt, output_path)
    logger.info('Successfully exported sites-only VCF')

    return output_path


if __name__ == '__main__':
    main()  # pylint: disable=E1120
