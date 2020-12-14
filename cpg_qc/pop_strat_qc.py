import logging
import pickle
from os.path import join
from typing import Optional, Tuple, List
import hail as hl
import pandas as pd

from gnomad.sample_qc.ancestry import run_pca_with_relateds, assign_population_pcs
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop
from gnomad.sample_qc.filtering import compute_qc_metrics_residuals, compute_stratified_metrics_filter, compute_stratified_sample_qc
from gnomad.sample_qc.pipeline import get_qc_mt

from gnomad.utils.file_utils import file_exists
from gnomad.utils.annotations import get_adj_expr

logger = logging.getLogger("cpg_qc_pca")


# Dictionary of file pointers that are used for checkpoints
TMP_PATH_MAP = {
    # Rows selected for PCA analysis
    'for_pca_mt':                     'for_pca.mt',

    'relatedness_ht':                 'relatedness.ht',
    'relatedness_pca_scores_ht':      'relatedness_pca_scores.ht',
    'pca_samples_rankings_ht':        'pca_samples_rankings.ht',
    'pop_pca_loadings_ht':            'pop_pca_loadings.ht',
    'pop_pca_eigenvalues_ht':         'pop_pca_eigenvalues.ht',
    'pop_pca_scores_ht':              'pop_pca_scores.ht',
    # Ranking of all release samples based on quality metrics. Used to remove relateds for release.
    'release_samples_ranking_ht':     'release_samples_ranking.ht',

    # Path to tab delimited file indicating inferred sample populations
    'pop_tsv':                        'RF_pop_assignments.txt.gz',
    # Path to RF model used for inferring sample populations
    'pop_rf':                         'pop.RF_fit.pickle',
}


def run_pop_strat_qc(
    mt: hl.MatrixTable,
    sample_df: pd.DataFrame,
    work_bucket: str,
    overwrite: bool,

    sex_ht: hl.Table,
    hard_filters_ht: hl.Table,
    sample_qc_ht: hl.Table,

    # Outputs:
    pca_related_samples_to_drop_ht_path: str,
    pop_ht_path: str,
    stratified_metrics_ht_path: str,
    regressed_metrics_ht_path: str,
    release_related_samples_to_drop_ht_path: str,

    # Parameters:
    kin_threshold: float,
    n_pcs: int,
    filtering_qc_metrics: List[str],
    min_pop_prob: float,
) -> Tuple[hl.Table, hl.Table, hl.Table, hl.Table, hl.Table]:

    tmp_path_map = {name: join(work_bucket, p) for name, p in TMP_PATH_MAP.items()}

    for_pca_mt = _make_mt_for_pca(
        mt,
        tmp_path_map['for_pca_mt'],
        overwrite
    )

    relatedness_ht = _compute_relatedness(
        for_pca_mt,
        sample_df,
        tmp_path_map['relatedness_pca_scores_ht'],
        tmp_path_map['relatedness_ht'],
        overwrite
    )

    pop_pca_scores_ht, pca_related_samples_to_drop_ht = _run_pca_ancestry_analysis(
        for_pca_mt=for_pca_mt,
        hard_filters_ht=hard_filters_ht,
        sex_ht=sex_ht,
        relatedness_ht=relatedness_ht,
        pca_samples_rankings_ht_path=tmp_path_map['pca_samples_rankings_ht'],
        pca_related_samples_to_drop_ht_path=pca_related_samples_to_drop_ht_path,
        pop_pca_scores_ht_path=tmp_path_map['pop_pca_scores_ht'],
        pop_pca_loadings_ht_path=tmp_path_map['pop_pca_loadings_ht'],
        pop_pca_eigenvalues_ht_path=tmp_path_map['pop_pca_eigenvalues_ht'],
        kin_threshold=kin_threshold,
        n_pcs=n_pcs,
        overwrite=overwrite,
    )

    pop_ht = _assign_pops(
        pop_pca_scores_ht,
        sample_df,
        pop_ht_path=pop_ht_path,
        pop_tsv_path=tmp_path_map['pop_tsv'],
        pop_rf_path=tmp_path_map['pop_rf'],
        min_prob=min_pop_prob,
        overwrite=overwrite,
    )

    stratified_metrics_ht = _compute_stratified_qc(
        sample_qc_ht,
        pop_ht,
        stratified_metrics_ht_path,
        filtering_qc_metrics=filtering_qc_metrics,
        overwrite=overwrite,
    )

    regressed_metrics_ht = _apply_regressed_filters(
        sample_qc_ht,
        pop_pca_scores_ht,
        regressed_metrics_ht_path,
        filtering_qc_metrics=filtering_qc_metrics,
        overwrite=overwrite
    )

    release_related_samples_to_drop_ht = _flag_related_samples(
        hard_filters_ht,
        sex_ht,
        relatedness_ht,
        regressed_metrics_ht,
        release_samples_rankings_ht_path=tmp_path_map['release_samples_ranking_ht'],
        release_related_samples_to_drop_ht_path=release_related_samples_to_drop_ht_path,
        kin_threshold=kin_threshold,
        overwrite=overwrite
    )

    return (
        pca_related_samples_to_drop_ht,
        pop_ht,
        stratified_metrics_ht,
        regressed_metrics_ht,
        release_related_samples_to_drop_ht
    )


def _make_mt_for_pca(
        mt: hl.MatrixTable,
        for_pca_mt_path: str,
        overwrite: bool
) -> hl.MatrixTable:
    """
    Create a new MatrixTable suitable for PCA analysis of relatedness
    and ancestry.
    * Strips unnesessary entry-level fields.
    * Kesps only bi-allelic SNPs.
    * Calls gnomad_methods's get_qc_mt() to filter rows further down based on:
      - the presence in problematic regions
      - callrate thresholds
      - MAF thresholds
      - inbreeding coefficient
      - allelic frequency thresholds
      - genotypes ADJ critetria (GQ>=20, DP>=10, AB>0.2 for hets)
    """
    if overwrite or not file_exists(for_pca_mt_path):
        logger.info('Making MatrixTable for PCA analysis')
        mt = mt.select_entries(
            'END',
            'LGT',
            GT=mt.LGT,
            adj=get_adj_expr(
                mt.LGT,
                mt.GQ,
                mt.DP,
                mt.LAD
            )
        )
        mt = mt.filter_rows(
            (hl.len(mt.alleles) == 2) &
            hl.is_snp(mt.alleles[0], mt.alleles[1])
        )
        mt = mt.naive_coalesce(5000)

        qc_mt = get_qc_mt(
            mt,
            adj_only=False,
            min_af=0.0,
            min_inbreeding_coeff_threshold=-0.025,
            min_hardy_weinberg_threshold=None,
            ld_r2=None,
            filter_lcr=False,
            filter_decoy=False,
            filter_segdup=False
        )
        qc_mt.write(for_pca_mt_path, overwrite=True)
    return hl.read_matrix_table(for_pca_mt_path)


def _compute_relatedness(
        for_pca_mt: hl.MatrixTable,
        sample_df: pd.DataFrame,
        relatedness_pca_scores_ht_path: str,
        relatedness_ht_path: str,
        overwrite: bool = False,
) -> str:

    if overwrite or not file_exists(relatedness_ht_path):
        logger.info('Running relatedness check')
        eig, scores, _ = hl.hwe_normalized_pca(
            for_pca_mt.GT,
            k=max(1, min(len(sample_df) // 3, 10)),
            compute_loadings=False
        )
        scores = scores.checkpoint(
            relatedness_pca_scores_ht_path,
            overwrite=overwrite, _read_if_exists=not overwrite)
        relatedness_ht = hl.pc_relate(
            for_pca_mt.GT,
            min_individual_maf=0.01,
            scores_expr=scores[for_pca_mt.col_key].scores,
            block_size=4096,
            min_kinship=0.05,
            statistics='all')

        # Converting keys for type struct{str} to str to align
        # with the rank_ht `s` key:
        relatedness_ht = relatedness_ht.key_by(
            i=relatedness_ht.i.s,
            j=relatedness_ht.j.s
        )

        relatedness_ht.write(relatedness_ht_path, overwrite=True)

    return hl.read_table(relatedness_ht_path)


def _run_pca_ancestry_analysis(
        for_pca_mt: hl.MatrixTable,
        hard_filters_ht: hl.Table,
        sex_ht: hl.Table,
        relatedness_ht: hl.Table,
        pca_samples_rankings_ht_path: str,
        pca_related_samples_to_drop_ht_path: str,
        pop_pca_scores_ht_path: str,
        pop_pca_loadings_ht_path: str,
        pop_pca_eigenvalues_ht_path: str,
        kin_threshold: float,
        n_pcs: int,
        overwrite: bool = False,
) -> Tuple[hl.Table, hl.Table]:

    if overwrite or not all(file_exists(fp) for fp in [
            pca_related_samples_to_drop_ht_path,
            pop_pca_scores_ht_path]):
        logger.info('Running PCA ancestry analysis')

        rank_ht = _compute_sample_rankings(
            hard_filters_ht,
            sex_ht,
            use_qc_metrics_filters=False,  # QC metrics filters ("regressed_metrics_ht")
                                           # do not exist at this point
        )
        rank_ht = rank_ht.checkpoint(
            pca_samples_rankings_ht_path,
            overwrite=overwrite, _read_if_exists=not overwrite)
        filtered_samples = hl.literal(rank_ht.aggregate(
            hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s))
        ))
        samples_to_drop = compute_related_samples_to_drop(
            relatedness_ht,
            rank_ht,
            kin_threshold,
            filtered_samples=filtered_samples
        )
        samples_to_drop.checkpoint(
            pca_related_samples_to_drop_ht_path,
            overwrite=overwrite, _read_if_exists=not overwrite)
        # Adjusting the number of principal components not to exceed the number of samples
        n_pcs = min(n_pcs, for_pca_mt.cols().count() - samples_to_drop.count())
        pop_pca_eignevalues, pop_pca_scores_ht, pop_pca_loadings_ht = run_pca_with_relateds(
            for_pca_mt,
            samples_to_drop,
            n_pcs=n_pcs
        )
        pop_pca_scores_ht.write(pop_pca_scores_ht_path, overwrite=True)
        pop_pca_loadings_ht.write(pop_pca_loadings_ht_path, overwrite=True)
        with hl.utils.hadoop_open(pop_pca_eigenvalues_ht_path, mode='w') as f:
            f.write(",".join([str(x) for x in pop_pca_eignevalues]))

    return (
        hl.read_table(pop_pca_scores_ht_path),
        hl.read_table(pca_related_samples_to_drop_ht_path),
    )


def _assign_pops(
        pop_pca_scores_ht: hl.Table,
        sample_df: pd.DataFrame,
        pop_ht_path: str,
        pop_tsv_path: str,
        pop_rf_path: str,
        min_prob: float,
        max_mislabeled_training_samples: int = 50,
        overwrite: bool = False,
) -> hl.Table:

    if overwrite or not all(file_exists(fp) for fp in [
            pop_ht_path, pop_tsv_path, pop_rf_path]):
        logger.info("Assigning global population labels")

        samples_with_pop_df = sample_df\
            [['sample', 'population']]\
            [pd.notna(sample_df['population'])]\
            .rename(columns={'sample': 's'})
        samples_with_pop_ht = hl.Table.from_pandas(samples_with_pop_df, key='s')
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=samples_with_pop_ht[pop_pca_scores_ht.key].population
        )

        def _run_assign_population_pcs(pop_pca_scores_ht, min_prob):
            examples_num = pop_pca_scores_ht.aggregate(
                hl.agg.count_where(hl.is_defined(pop_pca_scores_ht.training_pop))
            )
            logger.info(f'Running RF using {examples_num} training examples')
            pop_ht, pops_rf_model = assign_population_pcs(
                pop_pca_scores_ht,
                pc_cols=pop_pca_scores_ht.scores,
                known_col='training_pop',
                min_prob=min_prob
            )
            n_mislabeled_samples = pop_ht.aggregate(
                hl.agg.count_where(pop_ht.training_pop != pop_ht.pop))
            return pop_ht, pops_rf_model, n_mislabeled_samples

        pop_ht, pops_rf_model, n_mislabeled_samples = \
            _run_assign_population_pcs(pop_pca_scores_ht, min_prob)
        while n_mislabeled_samples > max_mislabeled_training_samples:
            logger.info(f"Found {n_mislabeled_samples} samples "
                        f"labeled differently from their known pop. "
                        f"Re-running without them.")

            pop_ht = pop_ht[pop_pca_scores_ht.key]
            pop_pca_scores_ht = pop_pca_scores_ht.annotate(
                training_pop=hl.or_missing(
                    (pop_ht.training_pop == pop_ht.pop),
                    pop_pca_scores_ht.training_pop
                )
            ).persist()

            pop_ht, pops_rf_model, n_mislabeled_samples = \
                _run_assign_population_pcs(pop_pca_scores_ht, min_prob)

        pop_ht = pop_ht.checkpoint(
            pop_ht_path,
            overwrite=overwrite, _read_if_exists=not overwrite)
        pc_cnt = min(hl.min(10, hl.len(pop_ht.pca_scores)).collect())
        pop_ht.transmute(
            **{f'PC{i + 1}': pop_ht.pca_scores[i] for i in range(pc_cnt)}
        ).export(pop_tsv_path)
        with hl.hadoop_open(pop_rf_path, 'wb') as out:
            pickle.dump(pops_rf_model, out)

    return hl.read_table(pop_ht_path)


def _compute_stratified_qc(
        sample_qc_ht: hl.Table,
        pop_ht: hl.Table,
        stratified_metrics_ht_path: str,
        filtering_qc_metrics: List[str],
        overwrite: bool = False,
) -> hl.Table:
    """
    Computes median, MAD, and upper and lower thresholds for each metric
    in `filtering_qc_metrics`, stratified by `pop` field in `pop_ht`
    """
    if overwrite or not file_exists(stratified_metrics_ht_path):
        logger.info(f'Computing stratified QC metrics filters using '
                    f'metrics: {", ".join(filtering_qc_metrics)}')

        sample_qc_ht = sample_qc_ht.annotate(
            qc_pop=pop_ht[sample_qc_ht.key].pop
        )
        stratified_metrics_ht = compute_stratified_metrics_filter(
            sample_qc_ht,
            qc_metrics={metric: sample_qc_ht.sample_qc[metric]
                        for metric in filtering_qc_metrics},
            strata={'qc_pop': sample_qc_ht.qc_pop},
            metric_threshold={'n_singleton': (4.0, 8.0)}
        )
        stratified_metrics_ht.write(stratified_metrics_ht_path,
                                    overwrite=True)
    return hl.read_table(stratified_metrics_ht_path)


def _flag_related_samples(
        hard_filters_ht: hl.Table,
        sex_ht: hl.Table,
        relatedness_ht: hl.Table,
        regressed_metrics_ht: hl.Table,
        release_samples_rankings_ht_path: str,
        release_related_samples_to_drop_ht_path: str,
        kin_threshold: float,
        overwrite: bool = False,
) -> hl.Table:

    if overwrite or not file_exists(release_related_samples_to_drop_ht_path):
        logger.info('Flagging related samples to drop')
        rank_ht = _compute_sample_rankings(
            hard_filters_ht,
            sex_ht,
            use_qc_metrics_filters=True,
            regressed_metrics_ht=regressed_metrics_ht,
        )
        rank_ht = rank_ht.checkpoint(
            release_samples_rankings_ht_path,
            overwrite=overwrite, _read_if_exists=not overwrite
        )
        filtered_samples = hl.literal(rank_ht.aggregate(
            hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s))))
        samples_to_drop = compute_related_samples_to_drop(
            relatedness_ht,
            rank_ht,
            kin_threshold=kin_threshold,
            filtered_samples=filtered_samples
        )
        samples_to_drop.write(release_related_samples_to_drop_ht_path, overwrite=True)

    return hl.read_table(release_related_samples_to_drop_ht_path)


def _compute_sample_rankings(
        hard_filtered_samples_ht: hl.Table,
        sex_ht: hl.Table,
        use_qc_metrics_filters: bool = False,
        regressed_metrics_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    :param hard_filtered_samples_ht: HT with hard_filters row field
    :param sex_ht: HT with chr20_mean_dp row field
    :param use_qc_metrics_filters: apply population-stratified QC filters
    :param regressed_metrics_ht: applied only if use_qc_metrics_filters
    :return: table ordered according to rank
    """
    ht = sex_ht.select(
        'chr20_mean_dp',
        filtered=hl.or_else(
            hl.len(hard_filtered_samples_ht[sex_ht.key].hard_filters) > 0,
            False)
    )
    if use_qc_metrics_filters:
        ht = ht.annotate(
            filtered=hl.cond(
                ht.filtered,
                True,
                hl.or_else(
                    hl.len(regressed_metrics_ht[ht.key].qc_metrics_filters) > 0,
                    False
                )
            )
        )

    ht = ht.order_by(
        ht.filtered,
        hl.desc(ht.chr20_mean_dp)
    ).add_index(name='rank')

    return ht.key_by('s').select('filtered', 'rank')


def _apply_regressed_filters(
        sample_qc_ht: hl.Table,
        ancestry_pca_scores_ht: hl.Table,
        regressed_metrics_ht_path: str,
        filtering_qc_metrics: List[str],
        overwrite: bool = False,
) -> hl.Table:
    """
    Compute QC metrics adjusted for popopulation
    """

    if overwrite or not file_exists(regressed_metrics_ht_path):
        logger.info('Compute QC metrics adjusted for popopulation')

        sample_qc_ht = sample_qc_ht.select(
            **sample_qc_ht.sample_qc,
            **ancestry_pca_scores_ht[sample_qc_ht.key],
            releasable=hl.bool(True)
        )
        residuals_ht = compute_qc_metrics_residuals(
            ht=sample_qc_ht,
            pc_scores=sample_qc_ht.scores,
            qc_metrics={metric: sample_qc_ht[metric] for metric in
                        filtering_qc_metrics},
            regression_sample_inclusion_expr=sample_qc_ht.releasable
        )

        stratified_metrics_ht = compute_stratified_metrics_filter(
            ht=residuals_ht,
            qc_metrics=dict(residuals_ht.row_value),
            metric_threshold={'n_singleton_residual': (4.0, 8.0)}
        )

        residuals_ht = residuals_ht.annotate(
            **stratified_metrics_ht[residuals_ht.key]
        )
        residuals_ht = residuals_ht.annotate_globals(
            **stratified_metrics_ht.index_globals()
        )

        residuals_ht.write(regressed_metrics_ht_path, overwrite=overwrite)

    return hl.read_table(regressed_metrics_ht_path)
