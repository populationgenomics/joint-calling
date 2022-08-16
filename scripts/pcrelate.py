#!/usr/bin/env python

"""
Run pcrelate.
"""

import logging

import click
import hail as hl

from large_cohort.query_utils import get_validation_callback, can_reuse, tmp_prefix

logger = logging.getLogger()
logger.setLevel('INFO')


@click.command()
@click.option(
    '--mt',
    'mt_path',
    required=True,
    callback=get_validation_callback(ext='mt', must_exist=True),
    help='path to a dense Matrix Table',
)
@click.option(
    '--out-ht',
    'out_ht_path',
    callback=get_validation_callback(ext='ht'),
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,missing-function-docstring
    mt_path: str,
    out_ht_path: str,
):
    """
    Writes table with the following structure:
    Row fields:
        'i': str
        'j': str
        'kin': float64
        'ibd0': float64
        'ibd1': float64
        'ibd2': float64
    Key: ['i', 'j']
    """
    if can_reuse(out_ht_path):
        return hl.read_table(out_ht_path)

    hl.init(default_reference='GRCh38')
    mt = hl.read_matrix_table(mt_path)
    mt = mt.select_entries('GT')

    logger.info('Running relatedness check')
    scores_ht_path = tmp_prefix / 'relatedness_pca_scores.ht'
    if can_reuse(scores_ht_path):
        scores_ht = hl.read_table(scores_ht_path)
    else:
        sample_num = mt.cols().count()
        _, scores_ht, _ = hl.hwe_normalized_pca(
            mt.GT, k=max(1, min(sample_num // 3, 10)), compute_loadings=False
        )
        scores_ht.checkpoint(scores_ht_path, overwrite=True)

    relatedness_ht = hl.pc_relate(
        mt.GT,
        min_individual_maf=0.01,
        scores_expr=scores_ht[mt.col_key].scores,
        block_size=4096,
        min_kinship=0.05,
    )

    # Converting keys for type struct{str} to str to align
    # with the rank_ht `s` key:
    relatedness_ht = relatedness_ht.key_by(i=relatedness_ht.i.s, j=relatedness_ht.j.s)
    return relatedness_ht.checkpoint(out_ht_path, overwrite=True)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
