#!/usr/bin/env python

"""
Convert somalier pairwise relatedness file into a Hail table similar to the 
one generated by hl.pc_relate, expected by the downstream analysis
"""

import logging
import pandas as pd
import click
import hail as hl

from joint_calling import utils
from joint_calling import _version


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--somalier-pairs-tsv',
    'somalier_pairs_tsv_path',
    required=True,
)
@click.option(
    '--out-relatedness-ht',
    'out_relatedness_ht_path',
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
    required=True,
    help='Hail billing account ID.',
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,missing-function-docstring
    somalier_pairs_tsv_path: str,
    out_relatedness_ht_path: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):
    utils.init_hail(__file__)

    if not utils.can_reuse(out_relatedness_ht_path, overwrite):
        pairs_df = pd.read_csv(somalier_pairs_tsv_path, delimiter='\t')
        pairs_df = pairs_df.rename(
            columns={'#sample_a': 'i', 'sample_b': 'j', 'relatedness': 'kin'}
        )
        ht = hl.Table.from_pandas(pairs_df).key_by('i', 'j')
        ht.write(out_relatedness_ht_path, overwrite=True)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
