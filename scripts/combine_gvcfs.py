#!/usr/bin/env python

"""
Combine a set of gVCFs and output a MatrixTable and a HailTable with metadata
"""

import os
from typing import List
import logging
import click
import hail as hl
from hail.experimental.vcf_combiner import vcf_combiner

from cpg_qc.utils import get_validation_callback, file_exists
from cpg_qc import utils
from cpg_qc import _version

logger = logging.getLogger('vcf_combiner')
logger.setLevel('INFO')

DEFAULT_REF = 'GRCh38'
# The target number of rows per partition during each round of merging
TARGET_RECORDS = 25_000


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--sample-map',
    'sample_map_csv_path',
    required=True,
    callback=get_validation_callback(ext='csv', must_exist=True),
    help='path to a per-sample data in a CSV file with '
    'a first line as a header. The only 2 required columns are `sample` '
    'and `gvcf`, in any order, possibly mixed with other columns.',
)
@click.option(
    '--out-mt',
    'out_mt_path',
    required=True,
    callback=get_validation_callback(ext='mt'),
    help='path to write the MatrixTable. Must have an .mt extention. '
    'Can be a Google Storage URL (i.e. start with `gs://`). '
    'An accompanying file with a `.metadata.ht` suffix will ne written '
    'at the same folder or bucket location, containing the same columns '
    'as the input sample map. This file is needed for further incremental '
    'extending of the matrix table using new GVCFs.',
)
@click.option(
    '--existing-mt',
    'existing_mt_path',
    callback=get_validation_callback(ext='mt', must_exist=True),
    help='optional path to an existing MatrixTable. Must have an .mt '
    'extention. Can be a Google Storage URL (i.e. start with `gs://`). '
    'If provided, will be read and used as a base to get extended with the '
    'samples in the input sample map. Can be read-only, as it will not '
    'be overwritten, instead the result will be written to the new location '
    'provided with --out-mt. An accompanying `.metadata.ht` file is expected '
    'to be present at the same folder or bucket location, containing the '
    'same set of samples, and the same columns as the input sample map',
)
@click.option(
    '--bucket',
    'work_bucket',
    required=True,
    help='path to folder for intermediate output. '
    'Can be a Google Storage URL (i.e. start with `gs://`).',
)
@click.option(
    '--local-tmp-dir',
    'local_tmp_dir',
    required=True,
    help='local directory for temporary files and Hail logs (must be local)',
)
@click.option(
    '--reuse',
    'reuse',
    is_flag=True,
    help='if an intermediate or a final file exists, reuse it instead of '
    'rerunning the code that generates it',
)
@click.option(
    '--hail-billing',
    'hail_billing',
    required=True,
    help='hail billing account ID',
)
def main(
    sample_map_csv_path: str,
    out_mt_path: str,
    existing_mt_path: str,
    work_bucket: str,
    local_tmp_dir: str,
    reuse: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):
    """
    Runs the Hail
    [vcf_combiner](https://hail.is/docs/0.2/experimental/vcf_combiner.html)
    using the GVCFs files specified in a `gvcf` column in the `sample_map_csv`
    CSV file as input, and generates a multi-sample Matrix Table in a sparse
    format, saved as `out_mt_path`. It also generates an accompanying table
    in an HT format with a `.metadata.ht` suffix, with the contents of the
    sample map, which can be used for incremental adding of new samples,
    as well as for running the QC.

    If `existing_mt_path` is provided, uses that matrix table as a base to
    extend with new samples. However, it will not overwrite `existing_mt_path`,
    and instead write the new table to `out_mt_path`. It would also combine
    the accompanying metadata HT tables and write the result with a
    `.metadata.ht` suffix.
    """
    utils.init_hail(
        name='combine_gvcfs',
        local_tmp_dir=local_tmp_dir,
    )

    new_metadata_ht = hl.import_table(sample_map_csv_path, delimiter=',', key='sample')

    if reuse and file_exists(existing_mt_path):
        logger.info(f'MatrixTable exists, reusing: {existing_mt_path}')
    else:
        logger.info(f'Combining new samples')
        new_mt_path = (
            os.path.join(work_bucket, 'new.mt') if existing_mt_path else out_mt_path
        )
        if reuse and file_exists(new_mt_path):
            logger.info(f'MatrixTable with new samples exists, reusing: {new_mt_path}')
        else:
            combine_gvcfs(
                gvcf_paths=new_metadata_ht.gvcf.collect(),
                out_mt_path=new_mt_path,
                work_bucket=work_bucket,
                overwrite=True,
            )
            logger.info(
                f'Written {new_metadata_ht.count()} new '
                f'samples into a MatrixTable {out_mt_path}'
            )
        if existing_mt_path:
            _combine_with_the_existing_mt(
                existing_mt=hl.read_matrix_table(existing_mt_path),
                new_mt_path=new_mt_path,
                out_mt_path=out_mt_path,
            )

    # Write metadata
    if existing_mt_path:
        existing_meta_ht_path = os.path.splitext(existing_mt_path)[0] + '.metadata.ht'
        existing_meta_ht = hl.read_table(existing_meta_ht_path)
        metadata_ht = existing_meta_ht.union(new_metadata_ht)
    else:
        metadata_ht = new_metadata_ht
    metadata_ht_path = os.path.splitext(out_mt_path)[0] + '.metadata.ht'
    if reuse and file_exists(metadata_ht_path):
        logger.info(f'Metadata table exists, reusing: {metadata_ht_path}')
    else:
        metadata_ht.write(metadata_ht_path, overwrite=True)
        logger.info(f'Written metadata table to {metadata_ht_path}')


def _combine_with_the_existing_mt(
    existing_mt: hl.MatrixTable,
    new_mt_path: str,  # passing as a path because we are going
    # to re-read it with different intervals
    out_mt_path: str,
):
    existing_mt = existing_mt.drop('gvcf_info')
    logger.info(
        f'Combining with the existing matrix table '
        f'({existing_mt.count_cols()} samples)'
    )
    intervals = vcf_combiner.calculate_new_intervals(
        hl.read_matrix_table(new_mt_path).rows(),
        n=TARGET_RECORDS,
        reference_genome=DEFAULT_REF,
    )
    new_mt = hl.read_matrix_table(new_mt_path, _intervals=intervals)
    new_mt = new_mt.drop('gvcf_info')
    out_mt = vcf_combiner.combine_gvcfs([existing_mt, new_mt])
    out_mt.write(out_mt_path, overwrite=True)


def combine_gvcfs(
    gvcf_paths: List[str], out_mt_path: str, work_bucket: str, overwrite: bool = True
):
    """
    Combine a set of GVCFs in one go
    """
    hl.experimental.run_combiner(
        gvcf_paths,
        out_file=out_mt_path,
        reference_genome=utils.DEFAULT_REF,
        use_genome_default_intervals=True,
        tmp_path=os.path.join(work_bucket, 'tmp'),
        overwrite=overwrite,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
