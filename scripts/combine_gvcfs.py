#!/usr/bin/env python

"""
Combine a set of GVCFs into a VDS.
"""

import logging

import click
import hail as hl
import pandas as pd
from cpg_utils import to_path
from cpg_utils.config import get_config
from hail.experimental.vcf_combiner.vcf_combiner import CombinerConfig
from hail.vds.combiner import new_combiner

from large_cohort.query_utils import check_duplicates, tmp_prefix

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
def main(
    cohort_tsv_path: str,
    out_vds_path: str,
):  # pylint: disable=missing-function-docstring
    hl.init(default_reference='GRCh38')

    branch_factor = (
        get_config()
        .get('combiner', {})
        .get('branch_factor', CombinerConfig.default_branch_factor)
    )
    batch_size = (
        get_config()
        .get('combiner', {})
        .get('batch_size', CombinerConfig.default_phase1_batch_size)
    )
    is_exome = get_config()['workflow']['sequencing_type'] == 'exome'

    logger.info(f'Combining GVCFs. Reading input paths from {cohort_tsv_path}')
    with to_path(cohort_tsv_path).open() as f:
        df = pd.read_table(f)
    sample_names = list(df.s)
    gvcf_paths = list(df.gvcf)
    check_duplicates(sample_names)
    check_duplicates(gvcf_paths)
    print(f'Combining {len(sample_names)} samples: {", ".join(sample_names)}')

    combiner = new_combiner(
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
    combiner.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
