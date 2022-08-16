#!/usr/bin/env python

"""
Select high-quality QC sites for PCA.
"""

import logging

import click
import hail as hl
from cpg_utils.hail_batch import reference_path
from hail import vds

from larcoh.query_utils import get_validation_callback

logger = logging.getLogger()
logger.setLevel('INFO')


@click.command()
@click.option(
    '--vds',
    'vds_path',
    required=True,
    callback=get_validation_callback(ext='vds', must_exist=True),
    help='path to the input dataset',
)
@click.option(
    '--out-mt',
    'out_mt_path',
    callback=get_validation_callback(ext='mt'),
    help='Path to resulting dense subset matrix table',
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,missing-function-docstring
    vds_path: str,
    out_mt_path: str,
):
    ds = vds.read_vds(vds_path)
    mt = vds.to_dense_mt(ds)
    sites_path = reference_path('ancestry_ht')
    sites_ht = hl.read_table(str(sites_path)).key_by('locus')

    mt = mt.filter_rows(
        (hl.len(mt.alleles) == 2)
        & hl.is_snp(mt.alleles[0], mt.alleles[1])
        & mt.locus.in_autosome()
        & hl.is_defined(sites_ht[mt.locus])
    )
    mt = mt.naive_coalesce(5000)
    mt.checkpoint(out_mt_path, overwrite=True)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
