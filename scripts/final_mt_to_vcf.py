#!/usr/bin/env python

"""
Convert matrix table to a VCF.
"""

import logging

import click
import hail as hl

from joint_calling.utils import get_validation_callback
from joint_calling import utils, _version

logger = logging.getLogger('joint-calling')
logger.setLevel('INFO')


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--mt',
    'mt_path',
    required=True,
    callback=get_validation_callback(ext='mt', must_exist=True),
    help='path to split matrix table',
)
@click.option(
    '--out-vcf',
    'out_vcf_path',
    required=True,
    help='path to write VCF',
)
@click.option(
    '--site-only',
    'site_only',
    is_flag=True,
    help='Do not export sample-level information. VCF will not have FORMAT and sample columns'
)
@click.option(
    '--local-tmp-dir',
    'local_tmp_dir',
    help='local directory for temporary files and Hail logs (must be local).',
)
@click.option(
    '--hail-billing',
    'hail_billing',
    help='Hail billing account ID.',
)
def main(
    mt_path: str,
    out_vcf_path: str,
    site_only: bool,
    local_tmp_dir: str,
    hail_billing: str,  # pylint: disable=unused-argument
):  # pylint: disable=missing-function-docstring

    utils.init_hail(__file__, local_tmp_dir)

    mt = hl.read_matrix_table(mt_path)
    write_vcf(mt, out_vcf_path, site_only=site_only)


def write_vcf(mt, out_vcf_path, site_only=False):
    """
    Convert final matrix table to a VCF.
    """
    mt = mt.drop('gvcf_info')

    if site_only:
        ds = mt.rows()
        annotate = ds.annotate
    else:
        ds = mt
        annotate = mt.annotate_rows
    
    ds = annotate(
        info=ds.info.annotate(
            InbreedingCoeff=ds.InbreedingCoeff,
            AC=ds.freq.AC, 
            AF=ds.freq.AF, 
            AN=ds.freq.AN, 
            homozygote_count=ds.freq.homozygote_count,
            faf95=ds.faf.faf95,
            faf99=ds.faf.faf99,
            popmax_AC=ds.popmax.AC, 
            popmax_AF=ds.popmax.AF, 
            popmax_AN=ds.popmax.AN, 
            popmax_homozygote_count=ds.popmax.homozygote_count, 
            popmax_pop=ds.popmax.pop,
            popmax_faf95=ds.popmax.faf95,
            nonsplit_alleles=ds.allele_data.nonsplit_alleles, 
            has_star=ds.allele_data.has_star, 
            variant_type=ds.allele_data.variant_type, 
            n_alt_alleles=ds.allele_data.n_alt_alleles, 
            allele_type=ds.allele_data.allele_type, 
            was_mixed=ds.allele_data.was_mixed, 
        )
    )
    hl.export_vcf(
        ds,
        out_vcf_path,
        tabix=True,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
