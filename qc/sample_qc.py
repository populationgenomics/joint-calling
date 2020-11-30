from typing import Optional

import hail as hl
from gnomad.resources.grch38 import (lcr_intervals, purcell_5k_intervals,
                                     telomeres_and_centromeres)
from gnomad.sample_qc.ancestry import (assign_population_pcs,
                                       run_pca_with_relateds)
from gnomad.sample_qc.filtering import (compute_qc_metrics_residuals,
                                        compute_stratified_metrics_filter,
                                        compute_stratified_sample_qc)
from gnomad.sample_qc.pipeline import annotate_sex, get_qc_mt
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop
from gnomad.sample_qc.sex import get_ploidy_cutoffs, get_sex_expr
from gnomad.utils.annotations import bi_allelic_expr, get_adj_expr
from gnomad.utils.filtering import add_filters_expr, filter_to_autosomes
from gnomad.utils.sparse_mt import densify_sites
from qc.resources import sex


def compute_sample_qc(
    mt: hl.MatrixTable,
    tmp_ht_prefix: Optional[str],
) -> hl.Table:
    mt = filter_to_autosomes(mt)
    mt = mt.filter_rows(~hl.is_defined(telomeres_and_centromeres.ht()[mt.locus]) & (hl.len(mt.alleles) > 1))
    mt = mt.select_entries('LGT')

    sample_qc_ht = compute_stratified_sample_qc(
        mt,
        strata={
            'bi_allelic': bi_allelic_expr(mt),
            'multi_allelic': ~bi_allelic_expr(mt)
        },
        tmp_ht_prefix=tmp_ht_prefix,
        gt_col='LGT'
    )

    # Remove annotations that cannot be computed from the sparse format
    sample_qc_ht = sample_qc_ht.annotate(
        **{
            x: sample_qc_ht[x].drop('n_called', 'n_not_called', 'call_rate')
            for x in sample_qc_ht.row_value
        }
    )
    return sample_qc_ht.repartition(100)


def compute_sex(mt: hl.MatrixTable,
                aaf_threshold=0.001, f_stat_cutoff=0.5) -> hl.Table:
    # Use AF from v3
    sex_ht = annotate_sex(
        mt,
        excluded_intervals=telomeres_and_centromeres.ht(),
        aaf_threshold=aaf_threshold,
        f_stat_cutoff=f_stat_cutoff,
        aaf_expr="AF",
        gt_expr="LGT",
    )

    return sex_ht













