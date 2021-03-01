"""
This script takes a path to a matrix table, and then:
 - exports a sites-only VCF using Hail Query,
 - run a pipeline on Cromwell
 - Collect the results

"""

import os
import logging

# import click

import hail as hl
import hailtop.batch as hb

from cpg_qc import _version
from cpg_qc.mt_to_vcf import mt_to_sites_only_mt


logger = logging.getLogger('vqsr_qc')
logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger.setLevel(logging.INFO)


# @click.command()
# @click.version_option(_version.__version__)
# TODO: add args here
def main_from_click(*args, **kwargs):
    """
    Driver function, but powered by @click
    """
    return main(*args, **kwargs)


def main(
    mt_path, output_dir: str, file_suffix: str, partitions: int, parallel: bool = None
):
    """
    Expects hail service to already be initialised
    """
    # do each step in hail batch

    logger.info(f'Loading matrix table from "{mt_path}"')
    mt = hl.read_matrix_table(mt_path)

    vcf_path = os.path.join(output_dir, f'{file_suffix}.sites.vcf.bgz')
    output_path = export_sites_only_vcf(
        mt=mt, output_path=vcf_path, partitions=partitions, parallel=parallel
    )

    # now run the WDL pipeline on Cromwell with inputs:
    inputs = {'vcf': output_path}
    print(inputs)

    return inputs


def export_sites_only_vcf(
    mt: hl.MatrixTable, output_path: str, partitions: int, parallel: bool
):
    """
    Take initial matrix table, convert to sites-only matrix table, then export to vcf
    """
    logger.info('Converting matrix table to sites-only matrix table')
    final_mt = mt_to_sites_only_mt(mt, partitions)

    # export vcf, and return the path

    logger.info(
        f"Exporting sites-only VCF to '{output_path}' to run in the VQSR pipeline"
    )
    hl.export_vcf(final_mt, output_path, parallel=parallel, tabix=True)
    logger.info('Successfully exported sites-only VCF')

    return output_path


if __name__ == '__main__':
    # main_from_click()

    hl.init(
        default_reference="GRCh38",
    )
    mt_path = "gs://cpg-fewgenomes-main/mt/50genomes.mt"
    output_dir = "gs://cpg-michael-dev-cromwell/vqsr-testing/"

    main(mt_path=mt_path, output_dir=output_dir, file_suffix="dev", partitions=1)
