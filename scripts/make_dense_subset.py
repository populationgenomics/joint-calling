#!/usr/bin/env python

"""
Combine with HGDP and select high-quality QC sites for PCA.

This script will produce two matrix tables:
1. `--out-hgdp-union-mt` - a union of all input dataset and HGDP-1KG samples,
   with rows corresponding to the rows shared with datasets, and fitered to
   high-quality ld-pruned sites best suited for PCA. The table is stipped off of
   all sample and entry annotations but GT, and with the HGDP-1KG sample annotations
   re-added as a `hgdp_1kg_metadata` column.
2. `--out-mt` - the input dataset with rows filtered down to the rows in 
   `--out-hgdp-union-mt` (so only high-quality sites shared with HGDP-1KG).
"""

import logging

import click
import hail as hl
from cpg_utils.hail_batch import reference_path
from hail import vds

from larcoh.query_utils import get_validation_callback

logger = logging.getLogger(__file__)


@click.command()
@click.option(
    '--vds',
    'vds_path',
    required=True,
    callback=get_validation_callback(ext='vds', must_exist=True),
    help='path to the input dataset',
)
@click.option(
    '--cohort-tsv',
    'cohort_tsv_path',
    callback=get_validation_callback(ext='tsv', must_exist=True),
    required=True,
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
    sites_path = reference_path('ancestry/v3-90k/pca_sites.ht')
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
