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
    '--n-pcs', 'n_pcs', type=int, help='number of PCs to compute for ancestry PCA.'
)
@click.option(
    '--pop',
    'pop',
    help='population label to subset the training dataset to',
)
@click.option(
    '--related-samples-to-drop-ht',
    'related_samples_to_drop_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
    help='Samples to remove from the analysis',
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
    related_samples_to_drop_ht_path: Optional[str],
    pop: Optional[str],
    out_eigenvalues_path: str,
    out_scores_ht_path: str,
    out_loadings_ht_path: str,
    tmp_bucket: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):
    utils.init_hail(__file__)

    if all(
        utils.can_reuse(fp, overwrite)
        for fp in [out_eigenvalues_path, out_scores_ht_path, out_loadings_ht_path]
    ):
        return

    mt = hl.read_matrix_table(mt_for_pca)

    if pop:
        # Get samples from the specified population only
        mt = mt.filter_cols(
            ~hl.is_defined(mt.hgdp_1kg_metadata)
            | (mt.hgdp_1kg_metadata.population_inference.pop == pop.lower())
        )

    related_samples_to_drop_ht = None
    if related_samples_to_drop_ht_path:
        related_samples_to_drop_ht = hl.read_table(related_samples_to_drop_ht_path)

    sqc.run_pca_ancestry_analysis(
        mt=mt,
        sample_to_drop_ht=related_samples_to_drop_ht,
        tmp_bucket=tmp_bucket,
        n_pcs=n_pcs,
        out_eigenvalues_path=out_eigenvalues_path,
        out_scores_ht_path=out_scores_ht_path,
        out_loadings_ht_path=out_loadings_ht_path,
        overwrite=overwrite,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
