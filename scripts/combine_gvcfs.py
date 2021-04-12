#!/usr/bin/env python

"""
Combine a set of gVCFs and output a MatrixTable and a HailTable with QC metadata
"""

import os
import subprocess
from typing import List
import logging
import click
import hail as hl
from hail.experimental.vcf_combiner import vcf_combiner

from joint_calling.utils import get_validation_callback, file_exists
from joint_calling import utils
from joint_calling import _version

logger = logging.getLogger('combine_gvcfs')
logger.setLevel('INFO')

DEFAULT_REF = 'GRCh38'
# The target number of rows per partition during each round of merging
TARGET_RECORDS = 25_000


@click.command()
@click.version_option(_version.__version__)
# @click.option(
#     '--dataset', 'dataset', required=True,
#     help='Dataset name'
# )
# @click.option(
#     '--test', 'test', required=True, is_flag=True,
#     help='Whether to use test or main bucket',
# )
# @click.option(
#     '--dataset-version', 'dataset_version', required=True,
#     help='Name or subfolder to find VCFs'
# )
@click.option('--bucket-with-vcfs', 'vcf_bucket')
@click.option('--qc-csv', 'qc_csv', required=True, help='File with QC metadata')
# @click.option(
#     '--sample-map',
#     'sample_map_csv_path',
#     required=True,
#     callback=get_validation_callback(ext='csv', must_exist=True),
#     help='path to a CSV file with per-sample data, where the '
#     'first line is a header. The only 2 required columns are `sample` '
#     '(the sample name) and `gvcf` (path to sample GVCF file) '
#     'in any order, possibly mixed with other columns.',
# )
@click.option(
    '--out-mt',
    'out_mt_path',
    required=True,
    callback=get_validation_callback(ext='mt'),
    help='path to write the MatrixTable. Must have an .mt extension. '
    'Can be a Google Storage URL (i.e. start with `gs://`). '
    'An accompanying file with a `.qc.ht` suffix will be written '
    'at the same folder or bucket location, containing the same columns '
    'as the input sample map. This file is needed for further incremental '
    'extending of the MatrixTable using new GVCFs.',
)
@click.option(
    '--existing-mt',
    'existing_mt_path',
    callback=get_validation_callback(ext='mt', must_exist=True),
    help='optional path to an existing MatrixTable. Must have an `.mt` '
    'extension. Can be a Google Storage URL (i.e. start with `gs://`). '
    'If provided, will be read and used as a base to get extended with the '
    'samples in the input sample map. Can be read-only, as it will not '
    'be overwritten, instead the result will be written to the new location '
    'provided with --out-mt. An accompanying `.qc.ht` file is expected '
    'to be present at the same folder or bucket location, containing the '
    'same set of samples, and the same columns as the input sample map.',
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
    help='local directory for temporary files and Hail logs (must be local).',
)
@click.option(
    '--reuse',
    'reuse',
    is_flag=True,
    help='if an intermediate or a final file exists, reuse it instead of '
    'rerunning the code that generates it.',
)
@click.option(
    '--hail-billing',
    'hail_billing',
    help='Hail billing account ID.',
)
def main(
    # sample_map_csv_path: str,
    # dataset: str,
    # test: bool,
    # dataset_version: str,
    vcf_bucket: str,
    qc_csv: str,
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
    using the GVCF files specified in a `gvcf` column in the `sample_map_csv`
    CSV file as input, and generates a multi-sample MatrixTable in a sparse
    format, saved as `out_mt_path`. It also generates an accompanying table
    in an HT format with a `.qc.ht` suffix, with the contents of the
    sample map, which can be used for incremental adding of new samples,
    as well as for running the QC.

    If `existing_mt_path` is provided, uses that MatrixTable as a base to
    extend with new samples. However, it will not overwrite `existing_mt_path`,
    and instead write the new table to `out_mt_path`. It would also combine
    the accompanying QC metadata HT tables and write the result with a
    `.qc.ht` suffix.
    """
    utils.init_hail('combine_gvcfs', local_tmp_dir)

    # sample.id,sample.sample_name,sample.flowcell_lane,sample.library_id,sample.platform,sample.centre,sample.reference_genome,raw_data.FREEMIX,raw_data.PlinkSex,raw_data.PCT_CHIMERAS,raw_data.PERCENT_DUPLICATION,raw_data.MEDIAN_INSERT_SIZE,raw_data.MEDIAN_COVERAGE
    # 613,TOB1529,ILLUMINA,HVTVGDSXY.1-2-3-4,LP9000039-NTP_H04,KCCG,hg38,0.0098939700,F(-1),0.023731,0.151555,412.0,31.0
    # 609,TOB1653,ILLUMINA,HVTVGDSXY.1-2-3-4,LP9000039-NTP_F03,KCCG,hg38,0.0060100100,F(-1),0.024802,0.165634,452.0,33.0
    # 604,TOB1764,ILLUMINA,HVTV7DSXY.1-2-3-4,LP9000037-NTP_B02,KCCG,hg38,0.0078874400,F(-1),0.01684,0.116911,413.0,43.0
    # 633,TOB1532,ILLUMINA,HVTVGDSXY.1-2-3-4,LP9000039-NTP_C05,KCCG,hg38,0.0121946000,F(-1),0.024425,0.151094,453.0,37.0
    new_gvcf_paths = [
        line.strip()
        for line in subprocess.check_output(
            f'gsutil ls \'{vcf_bucket}/*.g.vcf.gz\'', shell=True
        )
        .decode()
        .split()
    ]

    if reuse and file_exists(out_mt_path):
        logger.info(f'MatrixTable exists, reusing: {out_mt_path}')
    else:
        logger.info(f'Combining new samples')
        new_mt_path = (
            os.path.join(work_bucket, 'new.mt') if existing_mt_path else out_mt_path
        )
        combine_gvcfs(
            gvcf_paths=new_gvcf_paths,
            out_mt_path=new_mt_path,
            work_bucket=work_bucket,
            overwrite=True,
        )
        logger.info(f'Written samples into a MatrixTable {out_mt_path}')
        if existing_mt_path:
            _combine_with_the_existing_mt(
                existing_mt=hl.read_matrix_table(existing_mt_path),
                new_mt_path=new_mt_path,
                out_mt_path=out_mt_path,
            )

    # Write QC metadata
    qc_ht_path = os.path.splitext(out_mt_path)[0] + '.qc.ht'
    if reuse and file_exists(qc_ht_path):
        logger.info(f'QC table exists, reusing: {qc_ht_path}')
    else:
        new_qc_ht = hl.import_table(qc_csv, delimiter=',', impute=True)
        new_qc_ht = new_qc_ht.select(
            s=new_qc_ht['sample.sample_name'],
            freemix=new_qc_ht['raw_data.FREEMIX'],
            pct_chimeras=new_qc_ht['raw_data.PCT_CHIMERAS'],
            duplication=new_qc_ht['raw_data.PERCENT_DUPLICATION'],
            median_insert_size=new_qc_ht['raw_data.MEDIAN_INSERT_SIZE'],
            mean_coverage=new_qc_ht['raw_data.MEDIAN_COVERAGE'],
        ).key_by('s')

        if existing_mt_path:
            existing_qc_ht_path = os.path.splitext(existing_mt_path)[0] + '.qc.ht'
            existing_qc_ht = hl.read_table(existing_qc_ht_path)
            qc_ht = existing_qc_ht.union(new_qc_ht)
        else:
            qc_ht = new_qc_ht

        qc_ht.write(qc_ht_path, overwrite=True)
        logger.info(f'Written QC table to {qc_ht_path}')


def _combine_with_the_existing_mt(
    existing_mt: hl.MatrixTable,
    new_mt_path: str,  # passing as a path because we are going
    # to re-read it with different intervals
    out_mt_path: str,
):
    existing_mt = existing_mt.drop('gvcf_info')
    logger.info(
        f'Combining with the existing MatrixTable '
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
