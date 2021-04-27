"""Sample-level hard filtering based on Picard statistics"""

import logging
import hail as hl

from joint_calling.utils import file_exists


logger = logging.getLogger('sample_qc_hard_filtering')


def compute_hard_filters(
    mt: hl.MatrixTable,
    qc_ht: hl.MatrixTable,
    sex_ht: hl.Table,
    hail_sample_qc_ht: hl.Table,
    out_ht_path: str,
    cov_threshold: int,
    overwrite: bool = False,
) -> hl.Table:
    """
    Uses the sex imputation results, results of the sample_qc() run on
    bi-allelic variants, and Picard stats files specificed in `sample_df`,
    to apply filters to samples in `mt` and create a table with
    samples that fail at least one sampe.

    :param mt: input matrix table
    :param qc_ht: QC metadata generated by combine_gvcfs. Expected fields:
        contamination, alignment_summary_metrics, duplicate_metrics,
        insert_size_metrics, wgs_metrics (any of those are optional).
        Values must point to corresponding Picard stats files (see
        `_parse_metrics` for details)
    :param sex_ht: required fields: "sex_karyotype", "chr20_mean_dp"
    :param hail_sample_qc_ht: required fields:
        "bi_allelic_sample_qc { n_snp, n_singleton, r_het_hom_var }"
    :param out_ht_path: location to write the hard filtered samples Table
    :param cov_threshold: minimal chr20 coverage
    :param overwrite: overwrite checkpoints if they exist
    :return: table with samples failed the filters, and the following structure:
        's': str
        'hard_filters': set<str>  # a non-empty subset of { ambiguous_sex,
            sex_aneuploidy,  low_coverage, bad_biallelic_metrics, contamination,
            chimera, coverage, insert_size }
    """
    logger.info('Generating hard filters')
    if not overwrite and file_exists(out_ht_path):
        return hl.read_table(out_ht_path)

    ht = mt.cols()
    ht = ht.annotate(hard_filters=hl.empty_set(hl.tstr))

    # Helper function to add filters into the `hard_filters` set
    def add_filter(ht, expr, name):
        return ht.annotate(
            hard_filters=hl.if_else(
                expr & hl.is_defined(expr), ht.hard_filters.add(name), ht.hard_filters
            )
        )

    # Remove samples with ambiguous sex assignments
    ht = add_filter(ht, sex_ht[ht.key].sex_karyotype == 'ambiguous', 'ambiguous_sex')
    ht = add_filter(
        ht,
        ~hl.set({'ambiguous', 'XX', 'XY'}).contains(sex_ht[ht.key].sex_karyotype),
        'sex_aneuploidy',
    )

    # Remove low-coverage samples
    # chrom 20 coverage is computed to infer sex and used here
    ht = add_filter(ht, sex_ht[ht.key].chr20_mean_dp < cov_threshold, 'low_coverage')

    # Remove extreme raw bi-allelic sample QC outliers
    ht = add_filter(
        ht,
        (
            (hail_sample_qc_ht[ht.key].bi_allelic_sample_qc.n_snp > 8.0e6)
            | (hail_sample_qc_ht[ht.key].bi_allelic_sample_qc.n_snp < 2.4e6)
            | (hail_sample_qc_ht[ht.key].bi_allelic_sample_qc.n_singleton > 4e5)
            | (hail_sample_qc_ht[ht.key].bi_allelic_sample_qc.r_het_hom_var > 3.3)
        ),
        'bad_biallelic_metrics',
    )

    # Remove samples that fail picard metric thresholds, percents are not divided
    # by 100, e.g. 5% == 5.00, 5% != 0.05
    ht = add_filter(
        ht,
        hl.is_missing(qc_ht[ht.key].freemix) | (qc_ht[ht.key].freemix > 5.00),
        'contamination',
    )
    ht = add_filter(
        ht,
        hl.is_missing(qc_ht[ht.key].pct_chimeras) | (qc_ht[ht.key].pct_chimeras > 5.00),
        'chimera',
    )
    ht = add_filter(
        ht,
        hl.is_missing(qc_ht[ht.key].mean_coverage)
        | (qc_ht[ht.key].mean_coverage < 15.0),
        'coverage',
    )
    ht = add_filter(
        ht,
        hl.is_missing(qc_ht[ht.key].median_insert_size)
        | (qc_ht[ht.key].median_insert_size < 250.9),
        'insert_size',
    )
    ht = ht.filter(hl.len(ht.hard_filters) > 0)
    ht.write(out_ht_path, overwrite=True)
    return ht
