#!/usr/bin/env python

"""
Run sample QC on a MatrixTable, hard filter samples, add soft filter labels.
"""

import logging
import click

from joint_calling import utils
from joint_calling import sample_qc as sqc
from joint_calling import _version


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--mt',
    'mt_path',
    required=True,
    callback=utils.get_validation_callback(ext='mt', must_exist=True),
    help='path to the input Matrix Table. '
    'To generate it, run the `combine_gvcfs` script',
)
@click.option(
    '--meta-tsv',
    'meta_tsv_path',
    required=True,
    help='path to a CSV with QC metadata for the samples in the input Matrix Table. '
    'The following columns are expected: '
    's,freemix,pct_chimeras,duplication,median_insert_size. '
    'Must be keyed by "s".',
)
@click.option(
    '--filter-cutoffs-file',
    'filter_cutoffs_path',
    help=f'YAML file with filtering cutoffs',
)
@click.option(
    '--out-hard-filtered-samples-ht',
    'out_hard_filtered_samples_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=False),
    required=True,
)
@click.option(
    '--out-sex-ht',
    'out_sex_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=False),
    required=True,
)
@click.option(
    '--out-hail-sample-qc-ht',
    'out_hail_sample_qc_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=False),
    required=True,
)
@click.option(
    '--out-custom-qc-ht',
    'out_custom_qc_ht_path',
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
    '--local-tmp-dir',
    'local_tmp_dir',
    help='local directory for temporary files and Hail logs (must be local).',
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
    mt_path: str,
    meta_tsv_path: str,
    filter_cutoffs_path: str,
    out_hard_filtered_samples_ht_path: str,
    out_sex_ht_path: str,
    out_hail_sample_qc_ht_path: str,
    out_custom_qc_ht_path: str,
    tmp_bucket: str,
    local_tmp_dir: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):
    local_tmp_dir = utils.init_hail(__file__, local_tmp_dir)

    mt = utils.get_mt(mt_path, passing_sites_only=True)
    mt_split = utils.get_mt(mt_path, passing_sites_only=True, split=True)

    cutoffs_d = utils.get_filter_cutoffs(filter_cutoffs_path)

    input_metadata_ht = utils.parse_input_metadata(
        meta_tsv_path=meta_tsv_path,
        local_tmp_dir=local_tmp_dir,
    )

    # `hail_sample_qc_ht` row fields: sample_qc, bi_allelic_sample_qc
    hail_sample_qc_ht = sqc.compute_hail_sample_qc(
        mt=mt_split,
        tmp_bucket=tmp_bucket,
        out_ht_path=out_hail_sample_qc_ht_path,
        overwrite=overwrite,
    )

    # Calculate separately the number of non-gnomAD SNPs
    sqc.cpg_custom_metrics(
        split_mt=mt_split,
        tmp_bucket=tmp_bucket,
        out_ht_path=out_custom_qc_ht_path,
        overwrite=overwrite,
        count_gnomad_snps=True,
        count_chrx_het_hom=True,
    )

    # `sex_ht` row fields: is_female, chr20_mean_dp, sex_karyotype
    sex_ht = sqc.infer_sex(
        mt=mt,
        tmp_bucket=tmp_bucket,
        out_ht_path=out_sex_ht_path,
        overwrite=overwrite,
    )

    # `out_hard_filtered_samples_ht_path` row fields: hard_filters;
    # also includes only failed samples
    sqc.compute_hard_filters(
        mt=mt,
        picard_metrics_ht=input_metadata_ht,
        sex_ht=sex_ht,
        hail_sample_qc_ht=hail_sample_qc_ht,
        cutoffs_d=cutoffs_d['hardfiltering'],
        out_ht_path=out_hard_filtered_samples_ht_path,
        overwrite=overwrite,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
