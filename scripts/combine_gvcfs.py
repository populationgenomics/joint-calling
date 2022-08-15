#!/usr/bin/env python

"""
Combine a set of GVCFs into a MatrixTable
"""

import os
from typing import List, Optional
import logging
import shutil

import click
import pandas as pd

import hail as hl

from cpg_utils import to_path

from hail.vds.combiner import new_combiner
from hail.experimental.vcf_combiner.vcf_combiner import CombinerConfig

from larcoh.utils import check_duplicates


logger = logging.getLogger('combine_gvcfs')
logger.setLevel('INFO')


@click.command()
@click.option(
    '--cohort-tsv',
    'cohort_tsv_path',
    required=True,
    help='Path to a TSV files with sample metadata. Expected fields: s, gvcf',
)
@click.option(
    '--out-vds',
    'out_vds_path',
    required=True,
    help='path to write VDS',
)
@click.option(
    '--branch-factor',
    'branch_factor',
    required=False,
    type=click.INT,
    default=CombinerConfig.default_branch_factor,
    help='Combiner branch factor (gvcfs will be split into ~branch_factor batches)',
)
@click.option(
    '--batch-size',
    'batch_size',
    required=False,
    type=click.INT,
    default=CombinerConfig.default_phase1_batch_size,
    help='Combiner batch size for phase1 (size of each batch of gvcfs to process '
    'in 1 job)',
)
@click.option(
    '--tmp-prefix',
    'tmp_prefix',
    required=True,
    help='path to folder for intermediate output. '
    'Can be a Google Storage URL (i.e. start with `gs://`).',
)
@click.option(
    '--is-exome',
    'is_exome',
    is_flag=True,
    default=False,
)
def main(
    cohort_tsv_path: str,
    out_vds_path: str,
    branch_factor: int,
    batch_size: int,
    tmp_prefix: str,
    is_exome: bool,
):  # pylint: disable=missing-function-docstring
    hl.init(default_reference='GRCh38')
    logger.info(f'Combining GVCFs. Reading input paths from {cohort_tsv_path}')
    with to_path(cohort_tsv_path).open() as f:
        df = pd.read_table(f)
    sample_names = df.s
    gvcf_paths = df.gvcf
    check_duplicates(sample_names)
    check_duplicates(gvcf_paths)
    print(f'Combining {len(sample_names)} samples: {", ".join(sample_names)}')

    new_combiner(
        gvcf_paths=gvcf_paths,
        gvcf_sample_names=sample_names,
        # Header must be used with gvcf_sample_names, otherwise gvcf_sample_names
        # will be ignored. The first gvcf path works fine as a header because it will
        # be only read until the last line that begins with "#":
        gvcf_external_header=gvcf_paths[0],
        output_path=out_vds_path,
        reference_genome='GRCh38',
        use_genome_default_intervals=not is_exome,
        use_exome_default_intervals=is_exome,
        temp_path=tmp_prefix,
        force=True,
        branch_factor=branch_factor,
        batch_size=batch_size,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
