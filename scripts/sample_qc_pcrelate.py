#!/usr/bin/env python

"""
Run sample QC on a MatrixTable, hard filter samples, add soft filter labels.
"""

from os.path import join
import logging
from typing import Optional

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
    '--mt',
    'mt_path',
    required=True,
    callback=utils.get_validation_callback(ext='mt', must_exist=True),
    help='path to the Matrix Table ready for PC relate',
)
@click.option(
    '--out-relatedness-ht',
    'out_relatedness_ht_path',
)
@click.option(
    '--tmp-bucket',
    'tmp_bucket',
    required=True,
    help='path to write temporary intermediate output and checkpoints',
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
    mt_path: str,
    out_relatedness_ht_path: str,
    tmp_bucket: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):
    mt = utils.get_mt(mt_path, passing_sites_only=True, biallelic_snps_only=True)
    mt = mt.select_entries('END', GT=mt.LGT)

    _compute_relatedness(
        mt=mt,
        tmp_bucket=tmp_bucket,
        out_ht_path=out_relatedness_ht_path,
        overwrite=overwrite,
    )


def _compute_relatedness(
    mt: hl.MatrixTable,
    tmp_bucket: str,
    out_ht_path: Optional[str] = None,
    overwrite: bool = False,
) -> hl.Table:
    """
    :param mt: variants selected for PCA analysis (biallelic SNPs, annotated with GT)
    :param tmp_bucket: path to write checkpoints
    :param out_ht_path: path to write relatedness data
    :param overwrite: overwrite checkpoints if they exist
    :return: table with the following structure:
    Row fields:
        'i': str
        'j': str
        'kin': float64
        'ibd0': float64
        'ibd1': float64
        'ibd2': float64
    ----------------------------------------
    Key: ['i', 'j']
    """
    logger.info('Running relatedness check')
    out_ht_path = out_ht_path or join(tmp_bucket, 'relatedness.ht')
    if utils.can_reuse(out_ht_path, overwrite):
        return hl.read_table(out_ht_path)

    sample_num = mt.cols().count()

    _, scores, _ = hl.hwe_normalized_pca(
        mt.GT, k=max(1, min(sample_num // 3, 10)), compute_loadings=False
    )
    scores = scores.checkpoint(
        join(tmp_bucket, 'relatedness_pca_scores.ht'),
        overwrite=overwrite,
        _read_if_exists=not overwrite,
    )
    relatedness_ht = hl.pc_relate(
        mt.GT,
        min_individual_maf=0.01,
        scores_expr=scores[mt.col_key].scores,
        block_size=4096,
        min_kinship=0.05,
        statistics='all',
    )

    # Converting keys for type struct{str} to str to align
    # with the rank_ht `s` key:
    relatedness_ht = relatedness_ht.key_by(i=relatedness_ht.i.s, j=relatedness_ht.j.s)
    return relatedness_ht.checkpoint(out_ht_path, overwrite=True)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
