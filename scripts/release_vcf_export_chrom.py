import logging
import pickle

import click
import hail as hl

from gnomad.utils.vcf import (
    build_vcf_export_reference,
    rekey_new_reference,
)
from joint_calling import utils, _version
from joint_calling.utils import get_validation_callback


logger = logging.getLogger(__file__)
logger.setLevel('INFO')


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--ht',
    'ht_path',
    required=True,
    callback=get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--vcf-header-txt',
    'vcf_header_txt_path',
    required=True,
    callback=get_validation_callback(ext='txt', must_exist=True),
)
@click.option(
    '--out-vcf',
    'out_vcf_path',
    required=True,
    callback=get_validation_callback(ext='vcf.bgz'),
)
@click.option(
    '--name',
    'name',
    required=True,
    help='name of the dataset',
)
@click.option(
    '--chromosome',
    'chromosome',
    help='write a VCF for a chromosome',
)
@click.option(
    '--local-tmp-dir',
    'local_tmp_dir',
    help='local directory for temporary files and Hail logs (must be local).',
)
def main(
    ht_path: str,
    vcf_header_txt_path: str,
    out_vcf_path: str,
    name: str,
    chromosome: str,
    local_tmp_dir: str,
):  # pylint: disable=missing-function-docstring
    utils.init_hail(__file__, local_tmp_dir)

    ht = hl.read_table(ht_path)

    logger.info('Loading VCF header dict...')
    with hl.hadoop_open(vcf_header_txt_path, 'rb') as f:
        header_dict = pickle.load(f)

    if chromosome:
        logger.info(f'Exporting chromosome {chromosome}....')
        ht = hl.filter_intervals(ht, [hl.parse_locus_interval(chromosome)])

    export_reference = build_vcf_export_reference(name)

    hl.export_vcf(
        rekey_new_reference(ht, export_reference),
        out_vcf_path,
        metadata=header_dict,
        tabix=True,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
