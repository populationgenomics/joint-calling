#!/usr/bin/env python

"""
Write sample QC results in 2 formats: a sample-lebel Hail Table,
and a TSV file.
"""

from os.path import join, splitext
from typing import Optional
import logging
import click
import hail as hl

from gnomad.sample_qc.filtering import compute_stratified_sample_qc
from gnomad.sample_qc.pipeline import annotate_sex
from gnomad.utils.annotations import bi_allelic_expr
from gnomad.utils.filtering import (
    filter_to_autosomes,
    add_filters_expr,
)
from gnomad.utils.sparse_mt import filter_ref_blocks

from joint_calling.utils import file_exists
from joint_calling import sample_qc as sqc
from joint_calling import _version, utils


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--qc-bucket',
    'qc_bucket',
    help='Bucket with the results from sample_qc_calculate.py',
)
@click.option(
    '--age-csv', 'age_csv_path', help='CSV file with 2 columns: `sample` and `age`'
)
@click.option(
    '--hard-filtered-samples-ht',
    'hard_filtered_samples_ht_path',
)
@click.option(
    '--out-meta-ht',
    'out_meta_ht_path',
    required=True,
)
@click.option(
    '--out-meta-tsv',
    'out_meta_tsv_path',
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
    qc_bucket: str,
    age_csv_path: str,
    hard_filtered_samples_ht_path: str,
    out_meta_ht_path: str,
    out_meta_tsv_path: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):
    input_metadata_ht = hl.read_table(join(qc_bucket, sqc.INPUT_METADATA_HT_NAME))
    hail_sample_qc_ht = hl.read_table(join(qc_bucket, sqc.HAIL_SAMPLE_QC_HT_NAME))
    nongnomad_snps_ht = hl.read_table(join(qc_bucket, sqc.NONGNOMAD_SNPS_HT_NAME))
    sex_ht = hl.read_table(join(qc_bucket, sqc.SEX_HT_NAME))
    pop_ht = hl.read_table(join(qc_bucket, sqc.POP_HT_NAME))
    regressed_metrics_ht = hl.read_table(join(qc_bucket, sqc.REGRESSED_METRICS_HT_NAME))
    related_samples_to_drop_ht = hl.read_table(
        join(qc_bucket, sqc.RELATED_SAMPLES_TO_DROP_HT_NAME)
    )

    hard_filtered_samples_ht = hl.read_table(hard_filtered_samples_ht_path)

    if age_csv_path:
        age_ht = (
            hl.import_table(age_csv_path, delimiter=',', types={'age': 'float'})
            .rename({'TOBIID': 's'})
            .key_by('s')
        )
        sex_ht = sex_ht.annotate(
            age=age_ht[sex_ht.external_id].age,
        )

    # Combine all intermediate tables together
    _generate_metadata(
        sample_qc_ht=hail_sample_qc_ht,
        nongnomad_snps_ht=nongnomad_snps_ht,
        sex_ht=sex_ht,
        input_metadata_ht=input_metadata_ht,
        hard_filtered_samples_ht=hard_filtered_samples_ht,
        regressed_metrics_ht=regressed_metrics_ht,
        pop_ht=pop_ht,
        related_samples_to_drop_after_qc_ht=related_samples_to_drop_ht,
        out_ht_path=out_meta_ht_path,
        out_tsv_path=out_meta_tsv_path,
        overwrite=overwrite,
    )


def _compute_hail_sample_qc(
    mt: hl.MatrixTable,
    work_bucket: str,
    overwrite: bool = False,
) -> hl.Table:
    """
    Runs Hail hl.sample_qc() to generate a Hail Table with
    per-sample QC metrics.
    :param mt: input MatrixTable, multiallelics split, genotype in GT fields
    :param work_bucket: bucket path to write checkpoints
    :param overwrite: overwrite checkpoints if they exist
    :return: a Hail Table with the following row fields:
        'sample_qc': {
            <output of hl.sample_qc() (n_filtered, n_hom_ref, etc)>
        }
        'bi_allelic_sample_qc': {
            <output of hl.sample_qc() for bi-allelic variants only>
        }
        'multi_allelic_sample_qc': {
            <output of hl.sample_qc() for multi-allelic variants only>
        }
    """
    logger.info('Sample QC')
    out_ht_path = join(work_bucket, 'hail_sample_qc.ht')
    if not overwrite and file_exists(out_ht_path):
        logger.info(f'{out_ht_path} exists, reusing')
        return hl.read_table(out_ht_path)

    mt = filter_to_autosomes(mt)

    mt = mt.select_entries('GT')

    mt = filter_ref_blocks(mt)

    # Remove centromeres and telomeres incase they were included and any reference blocks
    tel_cent_ht = hl.read_table(utils.TEL_AND_CENT_HT_PATH)
    mt = mt.filter_rows(hl.is_missing(tel_cent_ht[mt.locus]))

    sample_qc_ht = compute_stratified_sample_qc(
        mt,
        strata={
            'bi_allelic': bi_allelic_expr(mt),
            'multi_allelic': ~bi_allelic_expr(mt),
        },
        tmp_ht_prefix=None,
    )

    # Remove annotations that cannot be computed from the sparse format
    sample_qc_ht = sample_qc_ht.annotate(
        **{
            x: sample_qc_ht[x].drop('n_called', 'n_not_called', 'call_rate')
            for x in sample_qc_ht.row_value
        }
    )
    sample_qc_ht = sample_qc_ht.repartition(100)
    sample_qc_ht.write(out_ht_path, overwrite=True)
    return sample_qc_ht


def _snps_not_in_gnomad(
    mt: hl.MatrixTable,
    work_bucket: str,
    overwrite: bool,
    gnomad_path: str = sqc.GNOMAD_HT_PATH,
) -> hl.Table:
    """
    Count the number of variants per sample that do not occur in gnomAD
    :param mt: MatrixTable, with multiallelics split, GT field for genotype
    :return: per-sample Table, with the only int field "nongnomad_snps"
    """
    # Filter to SNPs
    logger.info('Countings SNPs not in gnomAD')
    out_ht_path = join(work_bucket, 'notingnomad.ht')
    if not overwrite and file_exists(out_ht_path):
        return hl.read_table(out_ht_path)

    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

    # Get entries (table annotated with locus, allele, sample)
    ht = mt.entries()

    # Filter to those variants that are not in gnomad
    gnomad = hl.read_table(gnomad_path)
    ht = ht.key_by('locus', 'alleles').anti_join(gnomad)
    # Count non-gnomad variants for each sample
    stats_ht = ht.group_by(ht.s).aggregate(nongnomad_snps=hl.agg.count())
    stats_ht.write(out_ht_path, overwrite=True)
    return stats_ht


def _infer_sex(
    mt: hl.MatrixTable,
    work_bucket: str,
    overwrite: bool = False,
    target_regions: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Runs gnomad.sample_qc.annotate_sex() to infer sex
    :param mt: input MatrixTable
    :param work_bucket: bucket to write checkpoints
    :param overwrite: overwrite checkpoints if they exist
    :param target_regions: if exomes, it should correspond to the regions in the panel
    :return: a Hail Table with the following row fields:
        'is_female': bool
        'f_stat': float64
        'n_called': int64
        'expected_homs': float64
        'observed_homs': int64
        'chr20_mean_dp': float32
        'chrX_mean_dp': float32
        'chrY_mean_dp': float32
        'chrX_ploidy': float32
        'chrY_ploidy': float32
        'X_karyotype': str
        'Y_karyotype': str
        'sex_karyotype': str
    """
    logger.info('Inferring sex')
    out_ht_path = join(work_bucket, 'sex.ht')
    if not overwrite and file_exists(out_ht_path):
        return hl.read_table(out_ht_path)

    ht = annotate_sex(
        mt,
        excluded_intervals=hl.read_table(utils.TEL_AND_CENT_HT_PATH),
        included_intervals=target_regions,
        gt_expr='LGT',
    )

    ht.write(out_ht_path, overwrite=True)
    return ht


def _generate_metadata(
    sample_qc_ht: hl.Table,
    nongnomad_snps_ht: hl.Table,
    sex_ht: hl.Table,
    input_metadata_ht: hl.Table,
    hard_filtered_samples_ht: hl.Table,
    regressed_metrics_ht: hl.Table,
    pop_ht: hl.Table,
    related_samples_to_drop_after_qc_ht: hl.Table,
    out_ht_path: str,
    out_tsv_path: str,
    overwrite: bool = False,
) -> hl.Table:
    """
    Combine all intermediate tables into a single metadata table

    :param sample_qc_ht: table with a `bi_allelic_sample_qc` row field
    :param nongnomad_snps_ht: table with a `nongnomad_snps` row field
    :param sex_ht: table with the follwing row fields:
        `f_stat`, `n_called`, `expected_homs`, `observed_homs`
    :param input_metadata_ht: table with stats from the input metadata
        (includes Picard-tools stats and assigned population tags)
    :param hard_filtered_samples_ht: table with a `hard_filters` field
    :param regressed_metrics_ht: table with global fields
        `lms`, `qc_metrics_stats`;
        and metric row fields `n_snp_residual`, `r_ti_tv_residual`, etc;
        as well as corresponding `fail_n_snp_residual`, etc,
        and a `qc_metrics_filters` row field.
    :param pop_ht: table with the following row fields:
        `pop`, `prob_CEU`, `pca_scores`, `training_pop`.
    :param related_samples_to_drop_after_qc_ht:
        table with related samples, after re-ranking them based
        on population-stratified QC
    :param overwrite: overwrite checkpoints if they exist
    :return: table with relevant fields from input tables,
        annotated with the following row fields:
            'related_after_qc': bool = in `related_samples_to_drop_after_qc_ht`
            'related_before_qc': bool =
                in `related_samples_to_drop_before_qc_ht`
            'high_quality': bool =
                not `hard_filters_ht.hard_filters` &
                not `regressed_metrics_ht.qc_metrics_filters`
            'release': bool = high_quality & not `related`
    """
    logger.info('Generate the metadata HT and TSV files')

    if not overwrite and file_exists(out_ht_path):
        meta_ht = hl.read_table(out_ht_path)
    else:
        meta_ht = sex_ht.transmute(
            impute_sex_stats=hl.struct(
                f_stat=sex_ht.f_stat,
                n_called=sex_ht.n_called,
                expected_homs=sex_ht.expected_homs,
                observed_homs=sex_ht.observed_homs,
            )
        )

        meta_ht = meta_ht.annotate_globals(**regressed_metrics_ht.index_globals())

        meta_ht = meta_ht.annotate(
            sample_qc=sample_qc_ht[meta_ht.key].bi_allelic_sample_qc,
            nongnomad_snps=nongnomad_snps_ht[meta_ht.key].nongnomad_snps,
            **hard_filtered_samples_ht[meta_ht.key],
            **regressed_metrics_ht[meta_ht.key],
            **pop_ht[meta_ht.key],
            **input_metadata_ht[meta_ht.key],
            related=hl.is_defined(related_samples_to_drop_after_qc_ht[meta_ht.key]),
        )

        meta_ht = meta_ht.annotate(
            hard_filters=hl.or_else(meta_ht.hard_filters, hl.empty_set(hl.tstr)),
            release_filters=add_filters_expr(
                filters={'related': meta_ht.related},
                current_filters=meta_ht.hard_filters.union(meta_ht.qc_metrics_filters),
            ),
        )

        meta_ht = meta_ht.annotate(
            high_quality=(hl.len(meta_ht.hard_filters) == 0),
            release=hl.len(meta_ht.release_filters) == 0,
        )
        meta_ht.write(out_ht_path, overwrite=True)

    out_tsv_path = out_tsv_path or splitext(out_ht_path)[0] + '.tsv'
    if overwrite or not file_exists(out_tsv_path):
        n_pcs = meta_ht.aggregate(hl.agg.min(hl.len(meta_ht.pca_scores)))
        meta_ht = meta_ht.transmute(
            **{f'PC{i + 1}': meta_ht.pca_scores[i] for i in range(n_pcs)},
            hard_filters=hl.or_missing(
                hl.len(meta_ht.hard_filters) > 0, hl.delimit(meta_ht.hard_filters)
            ),
            qc_metrics_filters=hl.or_missing(
                hl.len(meta_ht.qc_metrics_filters) > 0,
                hl.delimit(meta_ht.qc_metrics_filters),
            ),
        )
        meta_ht.flatten().export(out_tsv_path)

    return meta_ht


if __name__ == '__main__':
    main()  # pylint: disable=E1120
