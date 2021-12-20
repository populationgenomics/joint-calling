#!/usr/bin/env python

"""
Combine with HGDP-1KG subset of gnomAD v3, and select high-quality QC sites for PCA.

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
from joint_calling import utils
from joint_calling import _version
from joint_calling.resources import ANCESTRY_SITES


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--mt',
    'mt_path',
    required=True,
    callback=utils.get_validation_callback(ext='mt', must_exist=True),
    help='path to the Matrix Table',
)
@click.option(
    '--meta-tsv',
    'meta_tsv_path',
    required=True,
    help='path to a CSV with QC and population metadata for the samples',
)
@click.option(
    '--pop',
    'pop',
)
@click.option(
    '--out-mt',
    'out_mt_path',
    callback=utils.get_validation_callback(ext='mt'),
    help='path to write the Matrix Table after subsetting to selected rows. '
    'The difference with --out-hgdp-union-mt is that it contains only the dataset '
    'samples',
)
@click.option(
    '--tmp-bucket',
    'tmp_bucket',
)
@click.option(
    '--overwrite/--reuse',
    'overwrite',
    is_flag=True,
    help='if an intermediate or a final file exists, skip running the code '
    'that generates it.',
)
@click.option(
    '--is-test',
    'is_test',
    is_flag=True,
    help='subset the gnomAD table to 20 samples',
)
@click.option(
    '--hail-billing',
    'hail_billing',
    required=True,
    help='Hail billing account ID.',
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,missing-function-docstring
    mt_path: str,
    meta_tsv_path: str,
    pop: Optional[str],
    out_mt_path: Optional[str],
    tmp_bucket: str,
    overwrite: bool,
    is_test: bool,  # pylint: disable=unused-argument
    hail_billing: str,  # pylint: disable=unused-argument
):
    local_tmp_dir = utils.init_hail(__file__)

    input_metadata_ht = utils.parse_input_metadata(meta_tsv_path, local_tmp_dir)

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
