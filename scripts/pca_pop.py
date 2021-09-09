#!/usr/bin/env python

"""
Runs PCA on one population specifically
"""

import logging
from os.path import join

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
    '--pca-hgdp-mt-path',
    'pca_with_hgdp_mt_path',
    required=True,
    callback=utils.get_validation_callback(ext='mt', must_exist=True),
    help='path to Matrix Table generated with sample_qc_subset_mt_for_pca.py',
)
@click.option(
    '--pop',
    'pop',
    default='nfe',
    help='population label to subset the training dataset to',
)
@click.option(
    '--n-pcs', 'n_pcs', default=20, help='number of PCs to compute for ancestry PCA.'
)
@click.option(
    '--meta-ht-path',
    'meta_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
    help='Sample-level metadata table with a `related` column field',
)
@click.option(
    '--out-analysis-bucket',
    'out_analysis_bucket',
    required=True,
    help='bucket location to write results to be used outside of the pipeline '
    '(usually a bucket with a read access for all users)',
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
    pca_with_hgdp_mt_path: str,
    pop: str,
    n_pcs: int,
    meta_ht_path: str,
    out_analysis_bucket: str,
    tmp_bucket: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):
    utils.init_hail(__file__)

    def _make_path(fname):
        return join(out_analysis_bucket, fname.replace('.ht', f'_{pop}.ht'))

    eigenvalues_ht_path = _make_path(sqc.PCA_EIGENVALUES_HT_NAME)
    scores_ht_path = _make_path(sqc.PCA_SCORES_HT_NAME)
    loadings_ht_path = _make_path(sqc.PCA_LOADINGS_HT_NAME)

    if all(
        utils.can_reuse(fp, overwrite)
        for fp in [eigenvalues_ht_path, scores_ht_path, loadings_ht_path]
    ):
        return

    mt = hl.read_matrix_table(pca_with_hgdp_mt_path)
    # Get samples from the specified population only
    mt = mt.filter_cols(
        ~hl.is_defined(mt.hgdp_1kg_metadata)
        | (mt.hgdp_1kg_metadata.population_inference.pop == pop.lower())
    )

    if meta_ht_path:
        meta_ht = hl.read_table(meta_ht_path)
        related_samples_to_drop_ht = meta_ht.filter(meta_ht.related)
    else:
        related_samples_to_drop_ht = None

    sqc.run_pca_ancestry_analysis(
        mt=mt,
        sample_to_drop_ht=related_samples_to_drop_ht,
        tmp_bucket=tmp_bucket,
        n_pcs=n_pcs,
        out_eigenvalues_ht_path=eigenvalues_ht_path,
        out_scores_ht_path=scores_ht_path,
        out_loadings_ht_path=loadings_ht_path,
        overwrite=overwrite,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
