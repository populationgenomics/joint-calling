import os
from collections import defaultdict
import pandas as pd
import hail as hl
import logging

from cpg_qc.utils import safe_mkdir, gs_cache_file, file_exists

logger = logging.getLogger("cpg_qc_hard_filtering")


def compute_hard_filters(
    mt: hl.MatrixTable,
    sex_ht: hl.Table,
    sample_qc_bi_allelic_ht: hl.Table,
    sample_df: pd.DataFrame,
    hard_filters_ht_path: str,
    work_bucket: str,
    local_tmp_dir: str,
    cov_threshold: int,
    overwrite: bool = False,
) -> hl.Table:

    if overwrite or not file_exists(hard_filters_ht_path):
        logger.info('Generating hard filters')

        # Parsing Picard metrics
        metrics_ht = _parse_metrics(sample_df, work_bucket, local_tmp_dir)

        ht = mt.cols()
        ht = ht.annotate(hard_filters=hl.empty_set(hl.tstr))

        def add_filter(ht, expr, name):
            return ht.annotate(hard_filters =
                hl.if_else(expr & hl.is_defined(expr),
                           ht.hard_filters.add(name),
                           ht.hard_filters)
            )

        # Remove samples with ambiguous sex assignments
        ht = add_filter(ht, sex_ht[ht.key].sex_karyotype == 'ambiguous',
                        "ambiguous_sex")
        ht = add_filter(ht, ~hl.set({'ambiguous', 'XX', 'XY'})
                        .contains(sex_ht[ht.key].sex_karyotype), "sex_aneuploidy")

        # Remove low-coverage samples
        # chrom 20 coverage is computed to infer sex and used here
        ht = add_filter(ht, sex_ht[ht.key].chr20_mean_dp < cov_threshold,
                        "low_coverage")

        # Remove extreme raw bi-allelic sample QC outliers
        ht = add_filter(ht, (
            (sample_qc_bi_allelic_ht[ht.key].sample_qc.n_snp > 3.75e6) |
            (sample_qc_bi_allelic_ht[ht.key].sample_qc.n_snp < 2.4e6) |
            (sample_qc_bi_allelic_ht[ht.key].sample_qc.n_singleton > 1e5) |
            (sample_qc_bi_allelic_ht[ht.key].sample_qc.r_het_hom_var > 3.3)
        ), "bad_qc_metrics")

        # Remove samples that fail picard metric thresholds, percents are not divided by 100, e.g. 5% == 5.00, %5 != 0.05
        ht = add_filter(ht, metrics_ht[ht.key].freemix > 5.00, "contamination")
        ht = add_filter(ht, metrics_ht[ht.key].pct_chimeras > 5.00, "chimera")
        ht = add_filter(ht, metrics_ht[ht.key].mean_coverage < 15, "coverage")
        ht = add_filter(ht, metrics_ht[ht.key].median_insert_size < 250,
                        "insert_size")
        ht = ht.filter(hl.len(ht.hard_filters) > 0)
        ht.write(hard_filters_ht_path, overwrite=True)

    return hl.read_table(hard_filters_ht_path)


def _parse_metrics(sample_df, work_bucket, local_tmp_dir):
    """
    * Contamination: freemix > 5% (`call-UnmappedBamToAlignedBam/UnmappedBamToAlignedBam/*/call-CheckContamination/*.selfSM`/`FREEMIX`)
    * Chimeras: > 5% (`call-AggregatedBamQC/AggregatedBamQC/*/call-CollectAggregationMetrics/*.alignment_summary_metrics`/`PCT_CHIMERAS`)
    * Duplication: > 30% (`call-UnmappedBamToAlignedBam/UnmappedBamToAlignedBam/*/call-MarkDuplicates/*.duplicate_metrics`/`PERCENT_DUPLICATION`)
    * Median insert size: < 250 (`call-AggregatedBamQC/AggregatedBamQC/*/call-CollectAggregationMetrics/*.insert_size_metrics`/`MEDIAN_INSERT_SIZE`)
    * Median coverage < 15X (`call-CollectWgsMetrics/*.wgs_metrics`/`MEDIAN_COVERAGE`)
    """
    data = defaultdict(list)

    for i, row in sample_df.iterrows():
        data['s'].append(row['sample'])

        contam = row.get('contamination')
        data['freemix'].append(
            _parse_picard_metric(contam, 'FREEMIX', local_tmp_dir))

        aln_sum_metrics = row.get('alignment_summary_metrics')
        data['pct_chimeras'].append(
            _parse_picard_metric(aln_sum_metrics, 'PCT_CHIMERAS', local_tmp_dir))

        dup_metrics = row.get('duplicate_metrics')
        data['duplication'].append(
            _parse_picard_metric(dup_metrics, 'PERCENT_DUPLICATION', local_tmp_dir))

        is_metrics = row.get('insert_size_metrics')
        data['median_insert_size'].append(
            _parse_picard_metric(is_metrics, 'MEDIAN_INSERT_SIZE', local_tmp_dir))

        wgs_metrics = row.get('wgs_metrics')
        data['mean_coverage'].append(
            _parse_picard_metric(wgs_metrics, 'MEDIAN_COVERAGE', local_tmp_dir))

    csv_path = os.path.join(work_bucket, 'sample_qc_metrics.tsv')
    pd.DataFrame.from_dict(data).to_csv(csv_path, sep='\t', index=False)
    ht = hl.import_table(
        csv_path,
        key='s',
        types={
            "s":                  hl.tstr,
            "freemix":            hl.tfloat32,
            "pct_chimeras":       hl.tfloat32,
            "duplication":        hl.tfloat32,
            "median_insert_size": hl.tint32,
            "mean_coverage":      hl.tint32,
        })
    return ht


def _parse_picard_metric(fpath, metric_name, local_tmp_dir):
    val = 'NA'
    if not fpath or pd.isnull(fpath):
        return val
    with open(gs_cache_file(fpath, local_tmp_dir)) as fh:
        idx = None
        for line in fh:
            if f"\t{metric_name}\t" in line:
                idx = line.split('\t').index(metric_name)
                continue
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
