#!/usr/bin/env python

"""
Generate a subset Matrix Table given a list of projects.

Will also remove reference blocks.
"""

import logging
from os.path import join
from typing import Optional, List

import click
import hail as hl

from joint_calling.utils import get_validation_callback
from joint_calling import utils, _version

logger = logging.getLogger('joint-calling')
logger.setLevel('INFO')


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--mt',
    'mt_path',
    required=True,
    callback=get_validation_callback(ext='mt', must_exist=True),
    help='path to final filtered split matrix table',
)
@click.option(
    '--subset-project',
    'subset_projects',
    multiple=True,
    required=True,
    help='Create a subset matrix table including only these projects'
)
@click.option(
    '--out-mt',
    'out_mt_path',
    required=True,
    callback=get_validation_callback(ext='mt', must_exist=False),
    help='path to write the subset matrix table',
)
@click.option(
    '--local-tmp-dir',
    'local_tmp_dir',
    help='local directory for temporary files and Hail logs (must be local).',
)
@click.option(
    '--hail-billing',
    'hail_billing',
    help='Hail billing account ID.',
)
def main(
    mt_path: str,
    subset_projects: List[str],
    out_mt_path: str,
    local_tmp_dir: str,
    hail_billing: str,  # pylint: disable=unused-argument
):  # pylint: disable=missing-function-docstring

    utils.init_hail(__file__, local_tmp_dir)

    # This will return a table with split multiallelics, and thus GT as a field
    mt = hl.read_matrix_table(mt_path)
    all_projects = list(set(mt.meta.project.collect()))
    logger.info(
        f'Full matrix table: {mt.count_cols()} samples, '
        f'{mt.count_rows() / 1_000_000}M rows, '
        f'{len(all_projects)} projects: {", ".join(all_projects)}'
    )

    subset_projects = list(set(subset_projects))
    diff_projects = set(subset_projects) - set(all_projects)
    if diff_projects:
        raise click.BadParameter(
            f'--subset-project values should be a subset of the existing matrix '
            f'table projects: mt.meta.project ({all_projects}). '
            f'The following projects are not in existing projects: {diff_projects} '
        )
    logger.info(f'Subsetting to {", ".join(subset_projects)}')

    mt = mt.filter_cols(hl.literal(subset_projects).contains(mt.meta.project))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    
    mt.write(out_mt_path, overwrite=True)
    mt = hl.read_matrix_table(out_mt_path)
    new_projects = list(set(mt.meta.project.collect()))
    logger.info(
        f'Subset matrix table: {mt.count_cols()} samples, {mt.count_rows()} rows, '
        f'{len(all_projects)} projects: {", ".join(new_projects)}'
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
