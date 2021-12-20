"""Functions to infer ancestry and perform ancestry-stratfied sample QC"""

import logging
import pickle
from os.path import join
from typing import Optional, List, Dict
import hail as hl
import pandas as pd

from gnomad.sample_qc.ancestry import run_pca_with_relateds, assign_population_pcs
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop
from gnomad.sample_qc.filtering import (
    compute_qc_metrics_residuals,
    compute_stratified_metrics_filter,
)
from gnomad.sample_qc.filtering import compute_stratified_sample_qc
from gnomad.sample_qc.pipeline import annotate_sex
from gnomad.utils.file_utils import file_exists
from gnomad.utils.annotations import bi_allelic_expr
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.sparse_mt import filter_ref_blocks

from joint_calling import utils, resources

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def compute_hail_sample_qc(
    mt: hl.MatrixTable,
    tmp_bucket: str,
    out_ht_path: Optional[str] = None,
    overwrite: bool = False,
) -> hl.Table:
    """
    Runs Hail hl.sample_qc() to generate a Hail Table with
    per-sample QC metrics.
    :param mt: input MatrixTable, multiallelics split, genotype in GT fields
    :param tmp_bucket: bucket path to write checkpoints
    :param out_ht_path: location to write the result to
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
    out_ht_path = out_ht_path or join(tmp_bucket, 'hail_sample_qc.ht')
    if utils.can_reuse(out_ht_path, overwrite):
        return hl.read_table(out_ht_path)

    mt = filter_to_autosomes(mt)

    mt = mt.select_entries('GT')

    mt = filter_ref_blocks(mt)

    # Remove centromeres and telomeres incase they were included and any reference blocks
    tel_cent_ht = hl.read_table(resources.TEL_AND_CENT_HT)
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
    return sample_qc_ht.checkpoint(out_ht_path, overwrite=True)


def cpg_custom_metrics(
    split_mt: hl.MatrixTable,
    tmp_bucket: str,
    out_ht_path: Optional[str] = None,
    overwrite: bool = False,
    gnomad_path: str = resources.GNOMAD_HT,
    count_gnomad_snps: bool = True,
    count_chrx_het_hom: bool = True,
) -> hl.Table:
    """
    Extra metrics for CPG reports.

    Resulting table fields:
        nongnomad_snps: int
        chrX_r_het_hom_var: float

    :param split_mt: MatrixTable, with multiallelics split, GT field for genotype
    :param tmp_bucket: bucket path to write checkpoints
    :param out_ht_path: location to write the result to
    :param overwrite: overwrite checkpoints if they exist
    :param gnomad_path: path to GnomAD Hail Table
    :return: per-sample Table.
    """
    out_ht_path = out_ht_path or join(tmp_bucket, 'custom_qc.ht')
    if utils.can_reuse(out_ht_path, overwrite):
        return hl.read_table(out_ht_path)

    def _count_gnomad_snps(mt, ht):
        logger.info(
            'Count the number of variants per sample that do not occur in gnomAD'
        )

        mt = mt.filter_rows(hl.len(mt.alleles) > 1)
        mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

        # Get entries (table annotated with locus, allele, sample)
        entries_ht = mt.entries()

        # Filter to those variants that are not in gnomad
        gnomad_ht = hl.read_table(gnomad_path)
        entries_ht = entries_ht.key_by('locus', 'alleles').anti_join(gnomad_ht)

        # Count non-gnomad variants for each sample
        cols_ht = entries_ht.group_by(entries_ht.s).aggregate(
            nongnomad_snps=hl.agg.count()
        )
        ht = ht.annotate(nongnomad_snps=cols_ht[ht.key].nongnomad_snps)
        return ht

    def _chrx_het_hom(mt, ht):
        logger.info('Counting het/hom ratio on chrX as an extra sex check')
        # Number of hom/hem on chrX, as an extra sex check measure
        chrx_mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX')])
        chrx_ht = hl.sample_qc(chrx_mt).cols()
        ht = ht.annotate(chrX_sample_qc=chrx_ht[ht.key].sample_qc)
        ht = ht.annotate(chrX_r_het_hom_var=ht.chrX_sample_qc.r_het_hom_var)
        return ht

    ht = split_mt.cols()
    if count_gnomad_snps:
        ht = _count_gnomad_snps(split_mt, ht)
    if count_chrx_het_hom:
        ht = _chrx_het_hom(split_mt, ht)
    return ht.checkpoint(out_ht_path, overwrite=True)


def infer_sex(
    mt: hl.MatrixTable,
    tmp_bucket: str,
    out_ht_path: Optional[str] = None,
    overwrite: bool = False,
    target_regions: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Runs gnomad.sample_qc.annotate_sex() to infer sex
    :param mt: input MatrixTable
    :param tmp_bucket: bucket to write checkpoints
    :param out_ht_path: location to write the result to
    :param overwrite: overwrite checkpoints if they exist
    :param target_regions: if exomes, it should correspond to the regions in the panel
    :return: a Table with the following row fields:
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
    out_ht_path = out_ht_path or join(tmp_bucket, 'sex.ht')
    if utils.can_reuse(out_ht_path, overwrite):
        return hl.read_table(out_ht_path)

    ht = annotate_sex(
        mt,
        excluded_intervals=hl.read_table(resources.TEL_AND_CENT_HT),
        included_intervals=target_regions,
        gt_expr='LGT',
    )
    return ht.checkpoint(out_ht_path, overwrite=True)


def run_pca_ancestry_analysis(
    mt: hl.MatrixTable,
    sample_to_drop_ht: Optional[hl.Table],
    n_pcs: int,
    out_eigenvalues_path: str,
    out_scores_ht_path: str,
    out_loadings_ht_path: str,
    overwrite: bool = False,
) -> hl.Table:
    """
    :param mt: variants usable for PCA analysis, combined with samples
        with known populations (HGDP, 1KG, etc)
    :param tmp_bucket: bucket path to write checkpoints
    :param sample_to_drop_ht: table with samples to drop based on
        previous relatedness analysis. With a `rank` row field
    :param n_pcs: maximum number of principal components
    :param overwrite: overwrite checkpoints if they exist
    :param out_eigenvalues_path: path to a txt file to write PCA eigenvalues
    :param out_scores_ht_path: path to write PCA scores
    :param out_loadings_ht_path: path to write PCA loadings
    :return: a Table with a row field:
        'scores': array<float64>
    """
    logger.info('Running PCA ancestry analysis')
    if all(
        not fp or utils.can_reuse(fp, overwrite)
        for fp in [
            out_eigenvalues_path,
            out_scores_ht_path,
            out_loadings_ht_path,
        ]
    ):
        return hl.read_table(out_scores_ht_path)

    # Adjusting the number of principal components not to exceed the
    # number of samples
    samples_to_drop_num = 0 if sample_to_drop_ht is None else sample_to_drop_ht.count()
    n_pcs = min(n_pcs, mt.cols().count() - samples_to_drop_num)
    eigenvalues, scores_ht, loadings_ht = run_pca_with_relateds(
        mt, sample_to_drop_ht, n_pcs=n_pcs
    )

    hl.Table.from_pandas(pd.DataFrame(eigenvalues)).export(out_eigenvalues_path)
    loadings_ht.write(out_loadings_ht_path, overwrite=True)
    scores_ht.write(out_scores_ht_path, overwrite=True)
    scores_ht = hl.read_table(out_scores_ht_path)
    return scores_ht


def infer_pop_labels(
    scores_ht: hl.Table,
    training_pop_ht: hl.Table,
    tmp_bucket: str,
    min_prob: float,
    max_mislabeled_training_samples: int = 50,
    n_pcs: int = 16,
    out_ht_path: Optional[str] = None,
    overwrite: bool = False,
) -> hl.Table:
    """
    Take population PCA results and training data, and run random forest
    to assign global population labels.

    :param scores_ht: output table of `_run_pca_ancestry_analysis()`
        with a row field 'scores': array<float64>
    :param training_pop_ht: table with samples with defined `training_pop`: str
    :param tmp_bucket: bucket to write checkpoints and intermediate files
    :param min_prob: min probability of belonging to a given population
        for the population to be set (otherwise set to `None`)
    :param max_mislabeled_training_samples: keep rerunning until the number
        of mislabeled samples is below this number
    :param n_pcs: Number of PCs to use in the RF
    :param out_ht_path: Path to write the resulting HT table
    :param overwrite: overwrite checkpoints if they exist
    :return: a Table with the following row fields, including `prob_<POP>`
        probabily fields for each population label:
        'training_pop': str
        'pca_scores': array<float64>
        'pop': str
        'prob_CEU': float64
        'prob_YRI': float64
        ... (prob_*: float64 for each population label)
    """
    logger.info('Assigning global population labels')
    out_ht_path = out_ht_path or join(tmp_bucket, 'inferred_pop.ht')
    if utils.can_reuse(out_ht_path, overwrite):
        return hl.read_table(out_ht_path)

    scores_ht = scores_ht.annotate(
        training_pop=training_pop_ht[scores_ht.key].training_pop
    )

    def _run_assign_population_pcs(pop_pca_scores_ht, min_prob):
        examples_num = pop_pca_scores_ht.aggregate(
            hl.agg.count_where(hl.is_defined(pop_pca_scores_ht.training_pop))
        )
        logger.info(f'Running RF using {examples_num} training examples')
        pop_ht, pops_rf_model = assign_population_pcs(
            pop_pca_scores_ht,
            pc_cols=pop_pca_scores_ht.scores[:n_pcs],
            known_col='training_pop',
            min_prob=min_prob,
        )
        n_mislabeled_samples = pop_ht.aggregate(
            hl.agg.count_where(pop_ht.training_pop != pop_ht.pop)
        )
        return pop_ht, pops_rf_model, n_mislabeled_samples

    pop_ht, pops_rf_model, n_mislabeled_samples = _run_assign_population_pcs(
        scores_ht, min_prob
    )
    while n_mislabeled_samples > max_mislabeled_training_samples:
        logger.info(
            f'Found {n_mislabeled_samples} samples '
            f'labeled differently from their known pop. '
            f'Re-running without them.'
        )

        pop_ht = pop_ht[scores_ht.key]
        pop_pca_scores_ht = scores_ht.annotate(
            training_pop=hl.or_missing(
                (pop_ht.training_pop == pop_ht.pop), scores_ht.training_pop
            )
        ).persist()

        pop_ht, pops_rf_model, n_mislabeled_samples = _run_assign_population_pcs(
            pop_pca_scores_ht, min_prob
        )

    # Writing a tab delimited file indicating inferred sample populations
    pop_tsv_file = join(tmp_bucket, 'RF_pop_assignments.txt.gz')
    if overwrite or not file_exists(pop_tsv_file):
        pc_cnt = min(hl.min(10, hl.len(pop_ht.pca_scores)).collect())
        pop_ht.transmute(
            **{f'PC{i + 1}': pop_ht.pca_scores[i] for i in range(pc_cnt)}
        ).export(pop_tsv_file)

    # Writing the RF model used for inferring sample populations
    pop_rf_file = join(tmp_bucket, 'pop.RF_fit.pickle')
    if overwrite or not file_exists(pop_rf_file):
        with hl.hadoop_open(pop_rf_file, 'wb') as out:
            pickle.dump(pops_rf_model, out)
    
    pop_ht = pop_ht.annotate(
        is_training=hl.is_defined(training_pop_ht[pop_ht.key])
    )
    return pop_ht.checkpoint(out_ht_path, overwrite=True)


def compute_stratified_qc(
    sample_qc_ht: hl.Table,
    pop_ht: hl.Table,
    work_bucket: str,
    filtering_qc_metrics: List[str],
    overwrite: bool = False,
) -> hl.Table:
    """
    Computes median, MAD, and upper and lower thresholds for each metric
    in `filtering_qc_metrics`, groupped by `pop` field in `pop_ht`

    :param sample_qc_ht: table with a row field
        `sample_qc` = struct { n_snp: int64, n_singleton: int64, ... }
    :param pop_ht: table with a `pop` row field
    :param work_bucket: bucket to write checkpoints
    :param filtering_qc_metrics: metrics to annotate with
    :param overwrite: overwrite checkpoints if they exist
    :return: table with the following structure:
        Global fields:
            'qc_metrics_stats': dict<tuple (
                str
                ), struct {
                    n_snp: struct {
                        median: int64,
                        mad: float64,
                        lower: float64,
                        upper: float64
                    },
                    n_singleton: struct {
                        ...
                    },
                    ...
            }
        Row fields:
            'fail_n_snp': bool
            'fail_n_singleton': bool
            ...
            'qc_metrics_filters': set<str>
    """
    logger.info(
        f'Computing stratified QC metrics filters using '
        f'metrics: {", ".join(filtering_qc_metrics)}'
    )

    sample_qc_ht = sample_qc_ht.annotate(qc_pop=pop_ht[sample_qc_ht.key].pop)
    stratified_metrics_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics={
            metric: sample_qc_ht.sample_qc[metric] for metric in filtering_qc_metrics
        },
        strata={'qc_pop': sample_qc_ht.qc_pop},
        metric_threshold={'n_singleton': (4.0, 8.0)},
    )
    return stratified_metrics_ht.checkpoint(
        join(work_bucket, 'stratified_metrics.ht'),
        overwrite=overwrite,
        _read_if_exists=not overwrite,
    )


def flag_related_samples(
    hard_filtered_samples_ht: hl.Table,
    sex_ht: hl.Table,
    relatedness_ht: hl.Table,
    regressed_metrics_ht: Optional[hl.Table],
    tmp_bucket: str,
    kin_threshold: float,
    out_ht_path: Optional[str] = None,
    overwrite: bool = False,
) -> hl.Table:
    """
    Flag samples to drop based on relatedness, so the final set
    has only unrelated samples, best quality one per family

    :param hard_filtered_samples_ht: table with failed samples
        and a `hard_filters` row field
    :param sex_ht: table with a `chr20_mean_dp` row field
    :param relatedness_ht: table keyed by exactly two fields (i and j)
        of the same type, identifying the pair of samples for each row
    :param regressed_metrics_ht: optional table with a `qc_metrics_filters`
        field calculated with _apply_regressed_filters() from PCA scores
    :param tmp_bucket: bucket to write checkpoints
    :param kin_threshold: kinship threshold to call two samples as related
    :param overwrite: overwrite checkpoints if they exist
    :param out_ht_path: path to write the resulting Table to
    :return: a table of the samples to drop along with their rank
        row field: 'rank': int64
    """
    label = 'final' if regressed_metrics_ht is not None else 'intermediate'
    logger.info(f'Flagging related samples to drop, {label}')
    out_ht_path = out_ht_path or join(tmp_bucket, f'{label}_related_samples_to_drop.ht')
    if utils.can_reuse(out_ht_path, overwrite):
        return hl.read_table(out_ht_path)

    rankings_ht_path = join(tmp_bucket, f'{label}_samples_rankings.ht')
    if utils.can_reuse(rankings_ht_path, overwrite):
        rank_ht = hl.read_table(rankings_ht_path)
    else:
        rank_ht = _compute_sample_rankings(
            hard_filtered_samples_ht,
            sex_ht,
            use_qc_metrics_filters=regressed_metrics_ht is not None,
            regressed_metrics_ht=regressed_metrics_ht,
        ).checkpoint(rankings_ht_path, overwrite=True)

    try:
        filtered_samples = hl.literal(
            rank_ht.aggregate(
                hl.agg.filter(rank_ht.filtered, hl.agg.collect(rank_ht.s))
            )
        )
    except hl.ExpressionException:
        # Hail doesn't handle it with `aggregate` when none of
        # the samples is 'filtered'
        filtered_samples = hl.empty_array('tstr')

    samples_to_drop_ht = compute_related_samples_to_drop(
        relatedness_ht,
        rank_ht,
        kin_threshold=kin_threshold,
        filtered_samples=filtered_samples,
    )
    return samples_to_drop_ht.checkpoint(out_ht_path, overwrite=True)


def _compute_sample_rankings(
    hard_filtered_samples_ht: hl.Table,
    sex_ht: hl.Table,
    use_qc_metrics_filters: bool = False,
    regressed_metrics_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Orders samples by hard filters and coverage and adds rank,
    which is the lower the better.

    :param hard_filtered_samples_ht: table with failed samples
        and a `hard_filters` row field
    :param sex_ht: table with a `chr20_mean_dp` row field
    :param use_qc_metrics_filters: apply population-stratified QC filters
    :param regressed_metrics_ht: table with a `qc_metrics_filters` field.
        Used only if `use_qc_metrics_filters` is True.
    :return: table ordered by rank, with the following row fields:
        `rank`, `filtered`
    """
    ht = sex_ht.drop(*list(sex_ht.globals.dtype.keys()))
    ht = ht.select(
        'chr20_mean_dp',
        filtered=hl.or_else(
            hl.len(hard_filtered_samples_ht[ht.key].hard_filters) > 0, False
        ),
    )
    if use_qc_metrics_filters and regressed_metrics_ht is not None:
        ht = ht.annotate(
            filtered=hl.cond(
                ht.filtered,
                True,
                hl.or_else(
                    hl.len(regressed_metrics_ht[ht.key].qc_metrics_filters) > 0, False
                ),
            )
        )

    ht = ht.order_by(ht.filtered, hl.desc(ht.chr20_mean_dp)).add_index(name='rank')
    return ht.key_by('s').select('filtered', 'rank')


def apply_regressed_filters(
    sample_qc_ht: hl.Table,
    pop_pca_scores_ht: hl.Table,
    tmp_bucket: str,
    out_ht_path: Optional[str] = None,
    overwrite: bool = False,
) -> hl.Table:
    """
    Re-compute QC metrics (with hl.sample_qc() - like n_snp, r_het_hom)
    per population, and adding "fail_*" row fields when a metric is below
    the the lower MAD threshold or higher the upper MAD threshold
    (see `compute_stratified_metrics_filter` for defaults)

    :param sample_qc_ht: table with a row field
       `bi_allelic_sample_qc` =
          struct { n_snp: int64, n_singleton: int64, ... }
    :param pop_pca_scores_ht: table with a `scores` row field
    :param tmp_bucket: bucket to write checkpoints
    :param out_ht_path: path to write the resulting table to
    :param overwrite: overwrite checkpoints if they exist
    :return: a table with the folliwing structure:
        Global fields:
            'lms': struct {
                n_snp: struct {
                    beta: array<float64>,
                    standard_error: array<float64>,
                    t_stat: array<float64>,
                    p_value: array<float64>,
                    multiple_standard_error: float64,
                    multiple_r_squared: float64,
                    adjusted_r_squared: float64,
                    f_stat: float64,
                    multiple_p_value: float64,
                    n: int32
                },
                n_singleton: struct {
                    ...
                },
                ...
            }
            'qc_metrics_stats': struct {
                n_snp_residual: struct {
                    median: float64,
                    mad: float64,
                    lower: float64,
                    upper: float64
                },
                n_singleton_residual: struct {
                    ...
                },
                ...
            }
        Row fields:
            's': str
            'n_snp_residual': float64
            'n_singleton_residual': float64
            ...
            'fail_n_snp_residual': bool
            'fail_n_singleton_residual': bool
            ...
            'qc_metrics_filters': set<str>
    """
    logger.info('Compute QC metrics adjusted for popopulation')
    out_ht_path = out_ht_path or join(tmp_bucket, 'regressed_metrics.ht')
    if utils.can_reuse(out_ht_path, overwrite):
        return hl.read_table(out_ht_path)

    sample_qc_ht = sample_qc_ht.select(
        **sample_qc_ht.bi_allelic_sample_qc,
        **pop_pca_scores_ht[sample_qc_ht.key],
        releasable=hl.bool(True),
    )

    filtering_qc_metrics = [
        'n_snp',
        'n_singleton',
        'r_ti_tv',
        'r_insertion_deletion',
        'n_insertion',
        'n_deletion',
        'r_het_hom_var',
        'n_het',
        'n_hom_var',
        'n_transition',
        'n_transversion',
    ]
    residuals_ht = compute_qc_metrics_residuals(
        ht=sample_qc_ht,
        pc_scores=sample_qc_ht.scores,
        qc_metrics={metric: sample_qc_ht[metric] for metric in filtering_qc_metrics},
        regression_sample_inclusion_expr=sample_qc_ht.releasable,
    )

    stratified_metrics_ht = compute_stratified_metrics_filter(
        ht=residuals_ht,
        qc_metrics=dict(residuals_ht.row_value),
        metric_threshold={'n_singleton_residual': (4.0, 8.0)},
    )

    residuals_ht = residuals_ht.annotate(**stratified_metrics_ht[residuals_ht.key])
    residuals_ht = residuals_ht.annotate_globals(
        **stratified_metrics_ht.index_globals()
    )

    return residuals_ht.checkpoint(out_ht_path, overwrite=True)


def compute_hard_filters(
    mt: hl.MatrixTable,
    picard_metrics_ht: hl.Table,
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
    :param picard_metrics_ht: table QC metrics from Picard tools. Expected
        fields: r_contamination, r_chimera, r_duplication, median_insert_size
    :param sex_ht: required fields: "sex_karyotype", "chr20_mean_dp"
    :param hail_sample_qc_ht: required fields:
        "bi_allelic_sample_qc { n_snp, n_singleton, r_het_hom_var }"
    :param out_ht_path: location to write the hard filtered samples Table
    :param cutoffs_d: a dictionary with hard-filtering thresholds
    :param overwrite: overwrite checkpoints if they exist
    :return: a table with samples failed the filters, and the following structure:
        's': str
        'hard_filters': set<str>  # a non-empty subset of { ambiguous_sex,
            sex_aneuploidy,  low_coverage, bad_biallelic_metrics, contamination,
            chimera, coverage, insert_size }
    """
    logger.info('Generating hard filters')
    if utils.can_reuse(out_ht_path, overwrite):
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
        hl.is_defined(picard_metrics_ht[ht.key].r_contamination)
        & (
            picard_metrics_ht[ht.key].r_contamination > cutoffs_d['max_r_contamination']
        ),
        'contamination',
    )
    ht = add_filter(
        ht,
        hl.is_defined(picard_metrics_ht[ht.key].r_chimera)
        & (picard_metrics_ht[ht.key].r_chimera > cutoffs_d['max_r_chimera']),
        'chimera',
    )
    ht = add_filter(
        ht,
        hl.is_defined(picard_metrics_ht[ht.key].r_duplication)
        & (picard_metrics_ht[ht.key].r_duplication > cutoffs_d['max_r_duplication']),
        'dup_rate',
    )
    ht = add_filter(
        ht,
        hl.is_defined(picard_metrics_ht[ht.key].median_insert_size)
        & (
            picard_metrics_ht[ht.key].median_insert_size
            < cutoffs_d['min_median_insert_size']
        ),
        'insert_size',
    )
    ht = ht.annotate_globals(hard_filter_cutoffs=hl.struct(**cutoffs_d))
    ht = ht.filter(hl.len(ht.hard_filters) > 0)
    return ht.checkpoint(out_ht_path, overwrite=True)
