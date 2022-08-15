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
from typing import Optional

import click
import hail as hl
from larcoh import utils

logger = logging.getLogger(__file__)


@click.command()
@click.option(
    '--mt',
    'mt_path',
    required=True,
    callback=utils.get_validation_callback(ext='mt', must_exist=True),
    help='path to the Matrix Table',
)
@click.option(
    '--cohort-tsv',
    'cohort_tsv_path',
    required=True,
    help='path to a CSV with QC metadata for the samples in the input Matrix Table. '
    'The following columns are expected: '
    's,freemix,pct_chimeras,duplication,insert_size. '
    'Must be keyed by "s".',
)
@click.option(
    '--out-mt',
    'out_mt_path',
    callback=utils.get_validation_callback(ext='mt'),
)
@click.option(
    '--out-combined-with-hgdp-mt',
    'out_combined_with_hgdp_mt_path',
    callback=utils.get_validation_callback(ext='mt'),
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,missing-function-docstring
    mt_path: str,
    cohort_tsv_path: str,
    out_mt_path: str,
    out_combined_with_hgdp_mt_path: str,
):
    utils.parse_input_metadata(meta_tsv_path, local_tmp_dir)

    sites_ht = hl.read_table(ANCESTRY_SITES).key_by('locus')

    mt = utils.get_mt(mt_path, passing_sites_only=True)
    mt = hl.experimental.densify(mt)
    mt = mt.select_entries(GT=mt.LGT).select_cols()
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
