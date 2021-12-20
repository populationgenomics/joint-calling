#!/usr/bin/env python

"""
Runs PCA, optionally on one population specifically
"""

import logging
from typing import Optional

import click
import hail as hl

from joint_calling import utils
from joint_calling import _version
from joint_calling import sample_qc as sqc


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--mt-for-pca',
    'mt_for_pca',
    required=True,
    callback=utils.get_validation_callback(ext='mt', must_exist=True),
    help='path to Matrix Table generated with sample_qc_subset_mt_for_pca.py',
)
@click.option(
    '--meta-tsv',
    'meta_tsv_path',
    callback=utils.get_validation_callback(ext='tsv', must_exist=True),
    required=True,
)
@click.option(
    '--n-pcs', 'n_pcs', type=int, help='number of PCs to compute for ancestry PCA.'
)
@click.option(
    '--related-samples-to-drop-ht',
    'related_samples_to_drop_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
    help='Samples to remove from the analysis',
)
@click.option(
    '--min-pop_prob', 
    'min_pop_prob',
    type=int,
    help='Minimal probability to infer population',
)
@click.option(
    '--subcontinental', 
    'subcontinental',
    is_flag=True,
    help='Infer subcontinental ancestry (default is continental level)',
)
@click.option(
    '--out-eigenvalues',
    'out_eigenvalues_path',
    callback=utils.get_validation_callback(ext='txt', must_exist=False),
)
@click.option(
    '--out-scores-ht',
    'out_scores_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=False),
)
@click.option(
    '--out-loadings-ht',
    'out_loadings_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=False),
)
@click.option(
    '--out-inferred-pop-ht',
    'out_inferred_pop_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=False),
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
    mt_for_pca: str,
    n_pcs: int,
    meta_tsv_path: str,
    related_samples_to_drop_ht_path: Optional[str],
    min_pop_prob: int,
    subcontinental: bool,
    out_eigenvalues_path: str,
    out_scores_ht_path: str,
    out_loadings_ht_path: str,
    out_inferred_pop_ht_path: str,
    tmp_bucket: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):
    local_tmp_dir = utils.init_hail(__file__)

    if all(utils.can_reuse(fp, overwrite) for fp in [
        out_eigenvalues_path, 
        out_scores_ht_path, 
        out_loadings_ht_path
    ]):
        return

    mt = hl.read_matrix_table(mt_for_pca)

    related_samples_to_drop_ht = None
    if related_samples_to_drop_ht_path:
        related_samples_to_drop_ht = hl.read_table(related_samples_to_drop_ht_path)

    scores_ht = sqc.run_pca_ancestry_analysis(
        mt=mt,
        sample_to_drop_ht=related_samples_to_drop_ht,
        n_pcs=n_pcs,
        out_eigenvalues_path=out_eigenvalues_path,
        out_scores_ht_path=out_scores_ht_path,
        out_loadings_ht_path=out_loadings_ht_path,
        overwrite=overwrite,
    )

    scores_ht = hl.read_table(scores_ht)
    
    meta_ht = utils.parse_input_metadata(meta_tsv_path, local_tmp_dir)

    tag = 'subcontinental_pop' if subcontinental else 'continental_pop'
    training_pop_ht = meta_ht.filter(
        hl.is_defined(meta_ht[tag]) & (meta_ht[tag] != '-')
    )
    training_pop_ht = training_pop_ht.annotate(
        training_pop=training_pop_ht[tag]
    )
    # Using calculated PCA scores as well as training samples with known
    # `population` tag, to assign population tags to remaining samples
    sqc.infer_pop_labels(
        scores_ht=scores_ht,
        training_pop_ht=training_pop_ht,
        tmp_bucket=tmp_bucket,
        min_prob=min_pop_prob,
        n_pcs=n_pcs,
        out_ht_path=out_inferred_pop_ht_path,
        overwrite=overwrite,
    )

# def _make_provided_pop_ht(
#     ht: hl.Table, 
#     subcontinental: bool = False
# ) -> hl.Table:
#     tag = 'subcontinental_pop' if subcontinental else 'continental_pop'
#     ht = ht.annotate(
#         pop=hl.if_else(ht[tag] != '-', ht[tag], ''),
#     )
#     return ht


if __name__ == '__main__':
    main()  # pylint: disable=E1120
