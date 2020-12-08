import logging
import os
import pickle
import traceback
from os.path import join
from typing import Optional, Tuple, List, Dict, Any
import hail as hl
import pandas as pd

from gnomad.sample_qc.ancestry import run_pca_with_relateds, assign_population_pcs
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop
from gnomad.sample_qc.filtering import compute_qc_metrics_residuals, compute_stratified_metrics_filter, compute_stratified_sample_qc
from gnomad.sample_qc.pipeline import annotate_sex, get_qc_mt

from gnomad.utils.file_utils import file_exists
from gnomad.utils.annotations import bi_allelic_expr, get_adj_expr

logger = logging.getLogger("cpg_qc")


TMP_PATH_MAP = {
    # Rows selected for PCA analysis
    'for_pca_mt':                     'for_pca.mt',

    'pca_scores_ht':                  'relatedness_pca_scores_ht',
    'relatedness_ht':                 'relatedness.ht',
    'relatedness_pca_scores_ht':      'relatedness_pca_scores.ht',
    'pca_samples_rankings_ht':        'pca_samples_rankings.ht',
    'ancestry_pca_loadings_ht':       'ancestry_pca_loadings.ht',
    'ancestry_pca_eigenvalues_ht':    'ancestry_pca_eigenvalues.ht',
    'ancestry_pca_scores_ht':         'ancestry_pca_scores_ht',
    # Ranking of all release samples based on quality metrics. Used to remove relateds for release.
    'release_samples_ranking_ht':     'release_samples_ranking.ht',

    # Path to tab delimited file indicating inferred sample populations
    'pop_tsv':                        'RF_pop_assignments.txt.gz',
    # Path to RF model used for inferring sample populations
    'pop_rf':                         'pop.RF_fit.pickle',
}


def pop_strat_qc(
        mt: hl.MatrixTable,
        sample_df: pd.DataFrame,
        work_bucket: str,
        overwrite: bool,

        sex_ht: hl.Table,
        hard_filters_ht: hl.Table,
        gnomad_meta_ht: Optional[hl.Table],
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
    ):

    path_map = {name: join(work_bucket, p)
                for name, p in TMP_PATH_MAP.items()}

    if overwrite or not file_exists(path_map['for_pca_mt']):
        """
        Creates a QC-ready MT by keeping:
        - Variants outside known problematic regions
        - Bi-allelic SNVs only
        - Variants passing hard thresholds
        - Variants passing the set call rate and MAF thresholds
        - Genotypes passing on gnomAD ADJ criteria (GQ>=20, DP>=10, AB>0.2 for hets)
        """
        logger.info('Making QC MT')
        for_pca_mt = compute_qc_mt(mt)
        for_pca_mt.write(path_map['for_pca_mt'], overwrite=True)

    if overwrite or not file_exists(path_map['relatedness_ht']):
        logger.info('Running relatedness')
        for_pca_mt = hl.read_matrix_table(path_map['for_pca_mt'])
        eig, scores, _ = hl.hwe_normalized_pca(
            for_pca_mt.GT,
            k=min(len(sample_df), 10),
            compute_loadings=False
        )
        scores = scores.checkpoint(path_map['relatedness_pca_scores_ht'],
                                   overwrite=overwrite, _read_if_exists=not overwrite)
        relatedness_ht = hl.pc_relate(
            for_pca_mt.GT,
            min_individual_maf=0.01,
            scores_expr=scores[for_pca_mt.col_key].scores,
            block_size=4096,
            min_kinship=0.05,
            statistics='all')
        try:
            relatedness_ht.write(path_map['relatedness_ht'], overwrite=True)
        except:
            # Can throw an error when the number of samples is small
            traceback.print_exc()

    if overwrite or not file_exists(path_map['ancestry_pca_scores_ht']):
        logger.info('Running PCA ancestry analysis')
        for_pca_mt = hl.read_matrix_table(path_map['for_pca_mt'])

        rank_ht = compute_sample_rankings(
            hard_filters_ht,
            sex_ht,
            use_qc_metrics_filters=False,  # QC metrics filters ("regressed_metrics_ht")
                                           # do not exist at this point
        )
        rank_ht = rank_ht.checkpoint(
            path_map['pca_samples_rankings_ht'],
            overwrite=overwrite, _read_if_exists=not overwrite)
        filtered_samples = hl.literal(rank_ht.aggregate(
            hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s))
        ))
        samples_to_drop = compute_related_samples_to_drop(
            hl.read_table(path_map['relatedness_ht']),
            rank_ht,
            kin_threshold,
            filtered_samples=filtered_samples
        )
        samples_to_drop.checkpoint(
            pca_related_samples_to_drop_ht_path,
            overwrite=overwrite, _read_if_exists=not overwrite)
        pop_pca_eignevalues, pop_pca_scores_ht, pop_pca_loadings_ht = run_pca_with_relateds(
            for_pca_mt,
            samples_to_drop,
            n_pcs=n_pcs
        )
        pop_pca_scores_ht.write(path_map['ancestry_pca_scores_ht'], overwrite=overwrite)
        pop_pca_loadings_ht.write(path_map['ancestry_pca_loadings_ht'],
                                  overwrite=overwrite)
        with hl.utils.hadoop_open(path_map['ancestry_pca_eigenvalues_ht'], mode='w') as f:
            f.write(",".join([str(x) for x in pop_pca_eignevalues]))

    if overwrite or not file_exists(pop_ht_path):
        pop_ht, pops_rf_model = assign_pops(
            hl.read_table(path_map['ancestry_pca_scores_ht']),
            gnomad_meta_ht,
            min_pop_prob
        )
        pop_ht = pop_ht.checkpoint(
            pop_ht_path,
            overwrite=overwrite, _read_if_exists=not overwrite)
        pop_ht.transmute(
            **{f'PC{i + 1}': pop_ht.pca_scores[i] for i in range(0, 10)}
        ).export(path_map['pop_tsv'])

        with hl.hadoop_open(path_map['pop_rf'], 'wb') as out:
            pickle.dump(pops_rf_model, out)

    if overwrite or not file_exists(stratified_metrics_ht_path):
        logger.info("Computing stratified QC metrics filters using "
                    "metrics: " + ", ".join(filtering_qc_metrics))
        pop_ht = hl.read_table(pop_ht_path)
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

    if overwrite or not file_exists(regressed_metrics_ht_path):
        # Computes qc_metrics adjusted for pop
        apply_regressed_filters(
            sample_qc_ht,
            hl.read_table(path_map['ancestry_pca_scores_ht']),
            gnomad_meta_ht,
            filtering_qc_metrics=filtering_qc_metrics,
        ).write(regressed_metrics_ht_path, overwrite=overwrite)

    if overwrite or not file_exists(release_related_samples_to_drop_ht_path):
        # Flags related samples to drop
        rank_ht = compute_sample_rankings(
            hard_filters_ht,
            sex_ht,
            use_qc_metrics_filters=True,
            regressed_metrics_ht=hl.read_table(path_map['regressed_metrics_ht']),
        )
        rank_ht = rank_ht.checkpoint(
            path_map['release_samples_rankings_ht'],
            overwrite=overwrite, _read_if_exists=not overwrite
        )
        filtered_samples = hl.literal(rank_ht.aggregate(
            hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s))))
        print(filtered_samples)
        samples_to_drop = compute_related_samples_to_drop(
            path_map['relatedness_ht'],
            rank_ht,
            kin_threshold,
            filtered_samples=filtered_samples
        )
        samples_to_drop.write(release_related_samples_to_drop_ht_path, overwrite=True)

def assign_pops(
        ancestry_pca_scores_ht: hl.Table,
        gnomad_meta_ht: Optional[hl.Table],
        min_prob: float,
        max_mislabeled_training_samples: int = 50
    ) -> Tuple[hl.Table, Any]:

    logger.info("Assigning global population labels")
    pop_pca_scores_ht = ancestry_pca_scores_ht.annotate(
        training_pop=(
            hl.case()
                .when(hl.is_defined(gnomad_meta_ht.project_pop), gnomad_meta_ht.project_pop)
                .when(gnomad_meta_ht.v2_pop != 'oth', gnomad_meta_ht.v2_pop)
                .or_missing()
        )
    )

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

    n_mislabeled_samples = pop_ht.aggregate(hl.agg.count_where(pop_ht.training_pop != pop_ht.pop))
    while n_mislabeled_samples > max_mislabeled_training_samples:
        logger.info(f"Found {n_mislabeled_samples} samples labeled differently from their known pop. Re-running without.")

        pop_ht = pop_ht[pop_pca_scores_ht.key]
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=hl.or_missing(
                (pop_ht.training_pop == pop_ht.pop),
                pop_pca_scores_ht.training_pop
            )
        ).persist()

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
        n_mislabeled_samples = pop_ht.aggregate(hl.agg.count_where(pop_ht.training_pop != pop_ht.pop))

    return pop_ht, pops_rf_model


def compute_qc_mt(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Create a new MatrixTable for QC purposes. Strips unnesessary
    entry-level fields, leaves only bi-allelic SNPs, and runs
    gnomad_methods's get_qc_mt() to filter rows further down
    based on callrate, inbreeding coeff, AF, MQ.
    The resulting set of rows is gonna be used for PC-relatedness.
    """
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
    return qc_mt


def compute_sample_rankings(
        hard_filtered_samples_ht: hl.Table,
        sex_ht: hl.Table,
        use_qc_metrics_filters: bool = False,
        regressed_metrics_ht: Optional[hl.Table] = None,  # applied only if use_qc_metrics_filters
) -> hl.Table:
    sex_ht = sex_ht.select(
        'chr20_mean_dp',
        filtered=hl.or_else(
            hl.len(hard_filtered_samples_ht[sex_ht.key].hard_filters) > 0,
            False)
    )
    if use_qc_metrics_filters:
        sex_ht = sex_ht.annotate(
            filtered=hl.cond(
                sex_ht.filtered,
                True,
                hl.or_else(
                    hl.len(regressed_metrics_ht[sex_ht.key].qc_metrics_filters) > 0,
                    False
                )
            )
        )

    sex_ht = sex_ht.order_by(
        sex_ht.filtered,
        hl.desc(sex_ht.chr20_mean_dp)
    ).add_index(name='rank')

    return sex_ht.key_by('s').select('filtered', 'rank')


def apply_regressed_filters(
        sample_qc_ht: hl.Table,
        ancestry_pca_scores_ht: hl.Table,
        project_meta_ht: hl.Table,
        filtering_qc_metrics: List[str],
) -> hl.Table:
    sample_qc_ht = sample_qc_ht.select(
        **sample_qc_ht.sample_qc,
        **ancestry_pca_scores_ht[sample_qc_ht.key],
        releasable=project_meta_ht[sample_qc_ht.key].releasable
    )
    residuals_ht = compute_qc_metrics_residuals(
        ht=sample_qc_ht,
        pc_scores=sample_qc_ht.scores,
        qc_metrics={metric: sample_qc_ht[metric] for metric in filtering_qc_metrics},
        regression_sample_inclusion_expr=sample_qc_ht.releasable
    )

    stratified_metrics_ht = compute_stratified_metrics_filter(
        ht=residuals_ht,
        qc_metrics=dict(residuals_ht.row_value),
        metric_threshold={'n_singleton_residual': (4.0, 8.0)}
    )

    residuals_ht = residuals_ht.annotate(
        **stratified_metrics_ht[sample_qc_ht.key]
    )
    residuals_ht = residuals_ht.annotate_globals(
        **stratified_metrics_ht.index_globals()
    )

    return residuals_ht
