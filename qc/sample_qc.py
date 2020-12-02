import os
from collections import defaultdict
from typing import Optional
import pandas as pd
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

from qc.utils import gs_download_file


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


def compute_hard_filters(mt: hl.MatrixTable,
                         cov_threshold: int,
                         biallelic_qc_ht: hl.Table,
                         sex_ht: hl.Table,
                         metrics_ht: hl.Table,
                         ) -> hl.Table:
    ht = mt.cols()
    hard_filters = dict()

    # Remove low-coverage samples
    cov_ht = sex_ht  # chrom 20 coverage is computed to infer sex and used here
    hard_filters['low_coverage'] = (cov_ht[ht.key].chr20_mean_dp < cov_threshold)

    # Remove extreme raw bi-allelic sample QC outliers
    hard_filters['bad_qc_metrics'] = (
            (biallelic_qc_ht.sample_qc.n_snp > 3.75e6) |
            (biallelic_qc_ht.sample_qc.n_snp < 2.4e6) |
            (biallelic_qc_ht.sample_qc.n_singleton > 1e5) |
            (biallelic_qc_ht.sample_qc.r_het_hom_var > 3.3)
    )

    # Remove samples with ambiguous sex assignments
    sex_ht = sex_ht[ht.key]
    hard_filters['ambiguous_sex'] = (sex_ht.sex_karyotype == 'Ambiguous')
    hard_filters['sex_aneuploidy'] = ~hl.set({'Ambiguous', 'XX', 'XY'}).contains(sex_ht.sex_karyotype)

    # Remove samples that fail picard metric thresholds, percents are not divided by 100, e.g. 5% == 5.00, %5 != 0.05
    hard_filters['contamination'] = metrics_ht.freemix > 5.00
    hard_filters['chimera'] = metrics_ht.pct_chimeras > 5.00
    hard_filters['coverage'] = metrics_ht.mean_coverage < 15
    hard_filters['insert_size'] = metrics_ht.median_insert_size < 250

    ht = ht.annotate(
        hard_filters=add_filters_expr(
            filters=hard_filters
        )
    )

    ht = ht.filter(hl.len(ht.hard_filters) > 0)
    return ht


def parse_metrics(sample_df, log_dirpath):
    """
	* Contamination: freemix > 5% (`call-UnmappedBamToAlignedBam/UnmappedBamToAlignedBam/*/call-CheckContamination/*.selfSM`/`FREEMIX`)
	* Chimeras: > 5% (`call-AggregatedBamQC/AggregatedBamQC/*/call-CollectAggregationMetrics/*.alignment_summary_metrics`/`PCT_CHIMERAS`)
	* Duplication: > 30% (`call-UnmappedBamToAlignedBam/UnmappedBamToAlignedBam/*/call-MarkDuplicates/*.duplicate_metrics`/`PERCENT_DUPLICATION`)
	* Median insert size: < 250 (`call-AggregatedBamQC/AggregatedBamQC/*/call-CollectAggregationMetrics/*.insert_size_metrics`/`MEDIAN_INSERT_SIZE`)
	* Median coverage < 15X (`call-CollectWgsMetrics/*.wgs_metrics`/`MEDIAN_COVERAGE`)
    """
    data = defaultdict(list)

    for i, row in sample_df.iterrows():
        sname = row['sample']

        contam = row.get('contamination')
        if not pd.isnull(contam):
            data['freemix'].append(_parse_picard_metric(contam, 'FREEMIX'))

        aln_sum_metrics = row.get('alignment_summary_metrics')
        if pd.isnull(aln_sum_metrics):
            data['chimera'].append(_parse_picard_metric(aln_sum_metrics, 'PCT_CHIMERAS'))

        dup_metrics = row.get('duplicate_metrics')
        if dup_metrics:
            data['duplication'].append(_parse_picard_metric(dup_metrics, 'PERCENT_DUPLICATION'))

        is_metrics = row.get('insert_size_metrics')
        if is_metrics:
            data['insert_size'].append(_parse_picard_metric(is_metrics, 'MEDIAN_INSERT_SIZE'))
        wgs_metrics = row.get('wgs_metrics')
        if wgs_metrics:
            data['coverage'].append(_parse_picard_metric(wgs_metrics, 'MEDIAN_COVERAGE'))

    csv_path = os.path.join(log_dirpath, 'metrics.tsv')
    pd.DataFrame.from_dict(data).to_csv(csv_path, sep='\t')
    ht = hl.import_table(csv_path, impute=True)
    return ht


def _parse_picard_metric(fpath, metric_name):
    val = None
    text = gs_download_file(fpath)
    idx = None
    for line in text.split('\n'):
        if f"\t{metric_name}\t" in line:
            idx = line.split('\t').index(metric_name)
        if idx is not None:
            val = line.split('\t')[idx]
            try:
                val = int(val)
            except ValueError:
                try:
                    val = float(val)
                except ValueError:
                    pass
                pass
            break
    return val








