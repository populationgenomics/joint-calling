#!/usr/bin/env python

"""
Run sample QC on a MatrixTable, hard filter samples, add soft filter labels.
"""

from os.path import join, basename
import subprocess
import logging
import click
import hail as hl
import pandas as pd

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
    '--meta-csv',
    'meta_csv_path',
    required=True,
    help='path to a CSV with QC and population metadata for the samples '
    'in the input Matrix Table. The following columns are expected: '
    's,population,gvcf,freemix,pct_chimeras,'
    'duplication,median_insert_size,mean_coverage. '
    'Must be keyed by "s". Samples with non-empty entries in '
    'the "population" column will be used to train the random forest '
    'for population inference of remaining samples. Other colums are '
    'used to apply QC hard filters to samples.',
)
@click.option(
    '--filter-cutoffs-file',
    'filter_cutoffs_path',
    help=f'YAML file with filtering cutoffs',
)
@click.option(
    '--relatedness-ht',
    'relatedness_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
    required=True,
)
@click.option(
    '--out-hard-filtered-samples-ht',
    'out_hard_filtered_samples_ht_path',
    callback=utils.get_validation_callback(ext='ht', must_exist=False),
    required=True,
)
@click.option(
    '--out-bucket', 'out_bucket', required=True, help='bucket location to write results'
)
@click.option(
    '--tmp-bucket',
    'tmp_bucket',
    required=True,
    help='path to write temporary intermediate output and checkpoints',
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
    '--n-pcs', 'n_pcs', default=30, help='number of PCs to compute for ancestry PCA.'
)
@click.option(
    '--target-bed',
    'target_bed',
    help='for exomes, target regions in a BED file format.',
)
@click.option(
    '--hail-billing',
    'hail_billing',
    required=True,
    help='Hail billing account ID.',
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,missing-function-docstring
    mt_path: str,
    meta_csv_path: str,
    filter_cutoffs_path: str,
    relatedness_ht_path: str,
    out_hard_filtered_samples_ht_path: str,
    out_bucket: str,
    tmp_bucket: str,
    local_tmp_dir: str,
    overwrite: bool,
    n_pcs: int,
    target_bed: str,
    hail_billing: str,  # pylint: disable=unused-argument
):
    input_metadata_ht_path = join(out_bucket, sqc.INPUT_METADATA_HT_NAME)
    hail_sample_qc_ht_path = join(out_bucket, sqc.HAIL_SAMPLE_QC_HT_NAME)
    nongnomad_snps_ht_path = join(out_bucket, sqc.NONGNOMAD_SNPS_HT_NAME)
    sex_ht_path = join(out_bucket, sqc.SEX_HT_NAME)
    pop_ht_path = join(out_bucket, sqc.POP_HT_NAME)
    regressed_metrics_ht_path = join(out_bucket, sqc.REGRESSED_METRICS_HT_NAME)

    local_tmp_dir = utils.init_hail('sample_qc', local_tmp_dir)

    mt = utils.get_mt(mt_path, passing_sites_only=True)
    mt_split = utils.get_mt(mt_path, passing_sites_only=True, split=True)
    mt_biall = utils.get_mt(mt_path, passing_sites_only=True, biallelic_snps_only=True)
    mt_biall = mt_biall.select_entries('END', GT=mt_biall.LGT)

    cutoffs_d = utils.get_filter_cutoffs(filter_cutoffs_path)

    input_metadata_ht = _parse_input_metadata(
        meta_csv_path=meta_csv_path,
        local_tmp_dir=local_tmp_dir,
        out_ht_path=input_metadata_ht_path,
    )

    # `hail_sample_qc_ht` row fields: sample_qc, bi_allelic_sample_qc
    hail_sample_qc_ht = sqc.compute_hail_sample_qc(
        mt=mt_split,
        tmp_bucket=tmp_bucket,
        out_ht_path=hail_sample_qc_ht_path,
        overwrite=overwrite,
    )

    # Calculate separately the number of non-gnomAD SNPs
    sqc.snps_not_in_gnomad(
        mt=mt_split,
        tmp_bucket=tmp_bucket,
        out_ht_path=nongnomad_snps_ht_path,
        overwrite=overwrite,
    )

    # `sex_ht` row fields: is_female, chr20_mean_dp, sex_karyotype
    sex_ht = sqc.infer_sex(
        mt=mt,
        tmp_bucket=tmp_bucket,
        out_ht_path=sex_ht_path,
        overwrite=overwrite,
        target_regions=hl.import_bed(target_bed) if target_bed else None,
    )

    # `hard_filtered_samples_ht` row fields: hard_filters;
    # also includes only failed samples
    hard_filtered_samples_ht = sqc.compute_hard_filters(
        mt=mt,
        picard_metrics_ht=input_metadata_ht,
        sex_ht=sex_ht,
        hail_sample_qc_ht=hail_sample_qc_ht,
        cutoffs_d=cutoffs_d['hardfiltering'],
        out_ht_path=out_hard_filtered_samples_ht_path,
        overwrite=overwrite,
    )

    # We don't want to include related samples into the
    # ancestry PCA analysis
    intermediate_related_samples_to_drop_ht = sqc.flag_related_samples(
        hard_filtered_samples_ht=hard_filtered_samples_ht,
        sex_ht=sex_ht,
        relatedness_ht=hl.read_table(relatedness_ht_path),
        regressed_metrics_ht=None,
        tmp_bucket=tmp_bucket,
        kin_threshold=cutoffs_d['pca']['max_kin'],
        overwrite=overwrite,
    )

    pop_pca_scores_ht = sqc.run_pca_ancestry_analysis(
        mt_biall=mt_biall,
        sample_to_drop_ht=intermediate_related_samples_to_drop_ht,
        tmp_bucket=tmp_bucket,
        n_pcs=n_pcs,
        overwrite=overwrite,
    )

    # Using calculated PCA scores as well as training samples with known
    # `population` tag, to assign population tags to remaining samples
    sqc.assign_pops(
        pop_pca_scores_ht=pop_pca_scores_ht,
        assigned_pop_ht=input_metadata_ht,
        tmp_bucket=tmp_bucket,
        min_prob=cutoffs_d['pca']['min_pop_prob'],
        n_pcs=n_pcs,
        out_ht_path=pop_ht_path,
        overwrite=overwrite,
    )

    # Re-computing QC metrics per population and annotating failing samples
    sqc.apply_regressed_filters(
        sample_qc_ht=hail_sample_qc_ht,
        pop_pca_scores_ht=pop_pca_scores_ht,
        tmp_bucket=tmp_bucket,
        out_ht_path=regressed_metrics_ht_path,
        overwrite=overwrite,
    )


def _parse_input_metadata(
    meta_csv_path: str,
    local_tmp_dir: str,
    out_ht_path: str,
) -> hl.Table:
    """
    Parse KCCG metadata (population and picard metrics)
    """
    local_csv_path = join(local_tmp_dir, basename(meta_csv_path))
    subprocess.run(
        f'gsutil cp {meta_csv_path} {local_csv_path}', check=False, shell=True
    )
    df = pd.read_table(local_csv_path)
    ht = hl.Table.from_pandas(df).key_by('s')
    return ht.checkpoint(out_ht_path, overwrite=True)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
