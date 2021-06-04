"""Sample-level hard filtering based on Picard statistics"""

import logging
from typing import Dict

import hail as hl

from joint_calling.utils import file_exists


logger = logging.getLogger('sample_qc_hard_filtering')


def compute_hard_filters(
    mt: hl.MatrixTable,
    input_meta_ht: hl.MatrixTable,
    sex_ht: hl.Table,
    hail_sample_qc_ht: hl.Table,
    cutoffs_d: Dict,
    out_ht_path: str,
    overwrite: bool = False,
) -> hl.Table:
    """
    Uses the sex imputation results, results of the sample_qc() run on
    bi-allelic variants, and stats by Picard tools specificed in `sample_df`,
    to apply filters to samples in `mt`, and create a table with
    only samples that fail at least one sample.

    :param mt: input matrix table
    :param input_meta_ht: QC metrics from Picard tools. Expected fields:
        r_contamination, r_chimera, r_duplication, median_insert_size
    :param sex_ht: required fields: "sex_karyotype", "chr20_mean_dp"
    :param hail_sample_qc_ht: required fields:
        "bi_allelic_sample_qc { n_snp, n_singleton, r_het_hom_var }"
    :param out_ht_path: location to write the hard filtered samples Table
    :param cutoffs_d: a dictionary with hard-filtering thresholds
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
    ht = add_filter(
        ht, sex_ht[ht.key].chr20_mean_dp < cutoffs_d['min_coverage'], 'low_coverage'
    )

    # Remove extreme raw bi-allelic sample QC outliers
    ht = add_filter(
        ht,
        (
            (
                hail_sample_qc_ht[ht.key].bi_allelic_sample_qc.n_snp
                > cutoffs_d['max_n_snps']
            )
            | (
                hail_sample_qc_ht[ht.key].bi_allelic_sample_qc.n_snp
                < cutoffs_d['min_n_snps']
            )
            | (
                hail_sample_qc_ht[ht.key].bi_allelic_sample_qc.n_singleton
                > cutoffs_d['max_n_singletons']
            )
            | (
                hail_sample_qc_ht[ht.key].bi_allelic_sample_qc.r_het_hom_var
                > cutoffs_d['max_r_het_hom']
            )
        ),
        'bad_biallelic_metrics',
    )

    ht = add_filter(
        ht,
        hl.is_missing(input_meta_ht[ht.key].r_contamination)
        | (input_meta_ht[ht.key].r_contamination > cutoffs_d['max_r_contamination']),
        'contamination',
    )
    ht = add_filter(
        ht,
        hl.is_missing(input_meta_ht[ht.key].r_chimera)
        | (input_meta_ht[ht.key].r_chimera > cutoffs_d['max_r_chimera']),
        'chimera',
    )
    ht = add_filter(
        ht,
        hl.is_missing(input_meta_ht[ht.key].r_duplication)
        | (input_meta_ht[ht.key].r_duplication > cutoffs_d['max_r_duplication']),
        'dup_rate',
    )
    ht = add_filter(
        ht,
        hl.is_missing(input_meta_ht[ht.key].median_insert_size)
        | (
            input_meta_ht[ht.key].median_insert_size
            < cutoffs_d['min_median_insert_size']
        ),
        'insert_size',
    )
    ht = ht.annotate_globals(hard_filter_cutoffs=hl.struct(**cutoffs_d))
    ht = ht.filter(hl.len(ht.hard_filters) > 0)
    ht.write(out_ht_path, overwrite=True)
    return ht
