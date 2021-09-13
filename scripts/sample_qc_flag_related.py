#!/usr/bin/env python

"""
Flag related samples to further exclude from the PCA analysis
"""

import logging
import click
import hail as hl

from joint_calling import sample_qc as sqc, utils
from joint_calling import _version


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--hard-filtered-samples-ht',
    'hard_filtered_samples_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--sex-ht',
    'sex_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
    help='Just for the `chr20_mean_dp` column',
)
@click.option(
    '--relatedness-ht',
    'relatedness_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--out-ht',
    'out_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=False),
)
@click.option(
    '--tmp-bucket',
    'tmp_bucket',
    required=True,
    help='path to write temporary intermediate output and checkpoints',
)
@click.option(
    '--max-kin',
    'max_kin',
    type=float,
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
    hard_filtered_samples_ht_path: str,
    sex_ht_path: str,
    relatedness_ht_path: str,
    out_ht_path: str,
    max_kin: int,
    tmp_bucket: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):
    if utils.can_reuse(out_ht_path, overwrite):
        return

    utils.init_hail(__file__)
    sqc.flag_related_samples(
        hard_filtered_samples_ht=hl.read_table(hard_filtered_samples_ht_path),
        sex_ht=hl.read_table(sex_ht_path),
        relatedness_ht=hl.read_table(relatedness_ht_path),
        regressed_metrics_ht=None,
        tmp_bucket=tmp_bucket,
        kin_threshold=max_kin,
        overwrite=overwrite,
        out_ht_path=out_ht_path,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
