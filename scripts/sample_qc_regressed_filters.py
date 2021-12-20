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
    '--hail-sample-qc-ht',
    'hail_sample_qc_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--filter-cutoffs-file',
    'filter_cutoffs_path',
    help=f'YAML file with filtering cutoffs',
)
@click.option(
    '--out-regressed-metrics-ht',
    'out_regressed_metrics_ht_path',
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
    hail_sample_qc_ht_path: str,
    out_regressed_metrics_ht_path: str,
    tmp_bucket: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):
    utils.init_hail(__file__)

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
