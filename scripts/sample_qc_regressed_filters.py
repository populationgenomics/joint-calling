#!/usr/bin/env python

"""
Run sample QC on a MatrixTable, hard filter samples, add soft filter labels.
"""

import logging
import click
import hail as hl

from joint_calling import utils
from joint_calling import sample_qc as sqc
from joint_calling import _version


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--pca-scores-ht',
    'pca_scores_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
    required=True,
)
@click.option(
    '--meta-tsv',
    'meta_tsv_path',
    callback=utils.get_validation_callback(ext='tsv', must_exist=True),
    required=True,
)
@click.option(
    '--hail-sample-qc-ht',
    'hail_sample_qc_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
    required=True,
)
@click.option(
    '--filter-cutoffs-file',
    'filter_cutoffs_path',
    help=f'YAML file with filtering cutoffs',
)
@click.option(
    '--n-pcs', 'n_pcs', default=20, help='number of PCs to compute for ancestry PCA.'
)
@click.option(
    '--out-regressed-metrics-ht',
    'out_regressed_metrics_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=False),
    required=True,
)
@click.option(
    '--out-inferred-pop-ht',
    'out_inferred_pop_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=False),
    required=True,
)
@click.option(
    '--tmp-bucket',
    'tmp_bucket',
    required=True,
    help='path to write temporary intermediate files and checkpoints '
    '(usually a temporary bucket)',
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
    pca_scores_ht_path: str,
    meta_tsv_path: str,
    hail_sample_qc_ht_path: str,
    filter_cutoffs_path: str,
    out_regressed_metrics_ht_path: str,
    out_inferred_pop_ht_path: str,
    tmp_bucket: str,
    overwrite: bool,
    n_pcs: int,
    hail_billing: str,  # pylint: disable=unused-argument
):
    local_tmp_dir = utils.init_hail(__file__)

    cutoffs_d = utils.get_filter_cutoffs(filter_cutoffs_path)
    
    meta_ht = utils.parse_input_metadata(meta_tsv_path, local_tmp_dir)

    # Using calculated PCA scores as well as training samples with known
    # `population` tag, to assign population tags to remaining samples
    sqc.infer_pop_labels(
        pop_pca_scores_ht=hl.read_table(pca_scores_ht_path),
        provided_pop_ht=meta_ht,
        tmp_bucket=tmp_bucket,
        min_prob=cutoffs_d['pca']['min_pop_prob'],
        n_pcs=n_pcs,
        out_ht_path=out_inferred_pop_ht_path,
        overwrite=overwrite,
    )

    # Re-computing QC metrics per population and annotating failing samples
    sqc.apply_regressed_filters(
        sample_qc_ht=hl.read_table(hail_sample_qc_ht_path),
        pop_pca_scores_ht=hl.read_table(pca_scores_ht_path),
        tmp_bucket=tmp_bucket,
        out_ht_path=out_regressed_metrics_ht_path,
        overwrite=overwrite,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
