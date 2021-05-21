#!/usr/bin/env python

"""
Run Random Forest variant QC
"""
from os.path import join
from typing import List, Optional, Tuple
import json
import logging
import uuid
import pyspark
import click
import hail as hl

from gnomad.resources.grch38.reference_data import (
    get_truth_ht,
    telomeres_and_centromeres,
)
from gnomad.variant_qc.pipeline import train_rf_model
from gnomad.variant_qc.random_forest import (
    apply_rf_model,
    get_rf_runs,
    get_run_data,
    load_model,
    median_impute_features,
    pretty_print_runs,
    save_model,
)

from joint_calling.utils import file_exists, get_validation_callback
from joint_calling import utils
from joint_calling import _version

logger = logging.getLogger('random_forest')
logger.setLevel('INFO')


FEATURES = [
    'allele_type',
    'AS_MQRankSum',
    'AS_pab_max',
    'AS_QD',
    'AS_ReadPosRankSum',
    'AS_SOR',
    'InbreedingCoeff',
    'n_alt_alleles',
    'variant_type',
]
INBREEDING_COEFF_HARD_CUTOFF = -0.3
INFO_FEATURES = [
    'AS_MQRankSum',
    'AS_pab_max',
    'AS_QD',
    'AS_ReadPosRankSum',
    'AS_SOR',
]
LABEL_COL = 'rf_label'
PREDICTION_COL = 'rf_prediction'
TRAIN_COL = 'rf_train'
TRUTH_DATA = ['hapmap', 'omni', 'mills', 'kgp_phase1_hc']


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--info-split-ht',
    'info_split_ht_path',
    required=True,
    callback=get_validation_callback(ext='ht', must_exist=True),
    help='path to info Table with split multiallelics '
    '(generated by generate_rf_annotations.py --split-multiallelic)',
)
@click.option(
    '--freq-ht',
    'freq_ht_path',
    required=True,
    callback=get_validation_callback(ext='ht', must_exist=True),
    help='path to a Table with InbreedingCoeff (generated by generate_freq_data.py)',
)
@click.option(
    '--fam-stats-ht',
    'fam_stats_ht_path',
    callback=get_validation_callback(ext='ht'),
    help='path to a Table with trio stats (generated by generate_qc_annotations.py)',
)
@click.option(
    '--allele-data-ht',
    'allele_data_ht_path',
    required=True,
    callback=get_validation_callback(ext='ht', must_exist=True),
    help='path to a Table with allele data (generated by generate_qc_annotations.py)',
)
@click.option(
    '--qc-ac-ht',
    'qc_ac_ht_path',
    required=True,
    callback=get_validation_callback(ext='ht', must_exist=True),
    help='path to a Table with allele counts (generated by generate_qc_annotations.py)',
)
@click.option(
    '--vqsr-filters-split-ht',
    'vqsr_filters_split_ht_path',
    callback=get_validation_callback(ext='ht'),
    help='Use VQSR training sites to train the RF (generated by load_vqsr.py)',
)
@click.option(
    '--out-annotations-ht',
    'out_annotations_ht_path',
    required=True,
    callback=get_validation_callback(ext='ht'),
    help='RF annotations',
)
@click.option(
    '--out-results-ht',
    'out_results_ht_path',
    required=True,
    callback=get_validation_callback(ext='ht'),
    help='RF results',
)
@click.option(
    '--out-model-id',
    'model_id',
    help='Model ID for the output models',
)
@click.option(
    '--bucket',
    'work_bucket',
    required=True,
    help='path to write intermediate output and checkpoints. '
    'Can be a Google Storage URL (i.e. start with `gs://`).',
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
    '--impute-features',
    'impute_features',
    is_flag=True,
    help='If set, feature imputation is performed',
)
@click.option(
    '--n-partitions',
    'n_partitions',
    type=click.INT,
    help='Desired number of partitions for annotated RF Table',
    default=5000,
)
@click.option(
    '--model-file',
    'model_file_path',
    help='Pre-created model file (requires --model-id and --model-training-ht)',
)
@click.option(
    '--model-training-ht',
    'model_training_ht_path',
    help='Pre-created model training file (requires --model-id and --model-file)',
)
# @click.option(
#     '--annotate-for-rf',
#     'do_annotate_for_rf',
#     is_flag=True,
#     help='Creates an annotated HT with features for RF.',
# )
# @click.option(
#     '--train-rf',
#     'do_train_rf',
#     is_flag=True,
#     help='Trains RF model.',
# )
# @click.option(
#     '--apply-rf',
#     'apply_rf_with_model_id',
#     help='Use this model id Applies RF model to the data.',
# )
# Random Forest parameters'
@click.option(
    '--fp-to-tp',
    'fp_to_tp',
    help='Ratio of FPs to TPs for training the RF model. If 0, '
    'all training examples are used. (default=1.0)',
    type=click.FLOAT,
    default=1.0,
)
@click.option(
    '--test-intervals',
    'test_intervals',
    help='The specified interval(s) will be held out for testing and '
    'evaluation only. (default to "chr20")',
    multiple=True,
    type=click.STRING,
    default=['chr20'],
)
@click.option(
    '--num-trees',
    'num_trees',
    type=click.INT,
    help='Number of trees in the RF model',
    default=500,
)
@click.option(
    '--max-depth',
    'max_depth',
    help='Maxmimum tree depth in the RF model',
    default=5,
    type=click.INT,
)
# Training data parameters
@click.option(
    '--use-adj-genotypes',
    'use_adj_genotypes',
    help='Use adj genotypes',
    is_flag=True,
)
# @click.option(
#     '--vqsr-model-id',
#     'vqsr_model_id',
#     help='If a VQSR model ID is provided the VQSR training annotations will be '
#          'used for training.',
#     default='vqsr_alleleSpecificTrans',
#     type=click.Choice(['vqsr_classic',
#                        'vqsr_alleleSpecific',
#                        'vqsr_alleleSpecificTrans'], case_sensitive=False),
# )
@click.option(
    '--no-transmitted-singletons',
    'no_transmitted_singletons',
    help='Do not use transmitted singletons for training.',
    is_flag=True,
)
@click.option(
    '--no-inbreeding-coeff',
    'no_inbreeding_coeff',
    help='Train RF without inbreeding coefficient as a feature.',
    is_flag=True,
)
@click.option(
    '--filter-centromere-telomere',
    'filter_centromere_telomere',
    help='Train RF without centromeres and telomeres.',
    is_flag=True,
)
def main(  # pylint: disable=too-many-arguments,too-many-locals
    info_split_ht_path: str,
    freq_ht_path: str,
    fam_stats_ht_path: str,
    allele_data_ht_path: str,
    qc_ac_ht_path: str,
    vqsr_filters_split_ht_path: Optional[str],
    out_results_ht_path: str,
    out_annotations_ht_path: str,
    model_id: str,
    work_bucket: str,
    local_tmp_dir: str,
    overwrite: bool,
    impute_features: bool,
    n_partitions: int,
    model_file_path: Optional[str],
    model_training_ht_path: Optional[str],
    fp_to_tp: Optional[float],
    test_intervals: Optional[List[str]],
    num_trees: Optional[int],
    max_depth: Optional[int],
    use_adj_genotypes: bool,
    # vqsr_model_id: str,
    no_transmitted_singletons: Optional[bool],
    no_inbreeding_coeff: Optional[bool],
    filter_centromere_telomere: Optional[bool],
):
    """
    Run variant QC on a MatrixTable, hard filter samples and add soft filter labels,
    output a sample-level Hail Table
    """
    local_tmp_dir = utils.init_hail('variant_qc_random_forest', local_tmp_dir)

    rf_json_path = join(work_bucket, 'rf_runs.json')
    rf_runs = get_rf_runs(rf_json_path)
    if rf_runs:
        logger.info(f'RF runs:')
        pretty_print_runs(rf_runs)

    if not overwrite and file_exists(out_annotations_ht_path):
        logger.info(f'Reusing {out_annotations_ht_path}')
        annotations_ht = hl.read_table(out_annotations_ht_path)
    else:
        annotations_ht = create_rf_ht(
            info_ht=hl.read_table(info_split_ht_path),
            inbreeding_ht=hl.read_table(freq_ht_path),
            trio_stats_ht=hl.read_table(fam_stats_ht_path)
            if fam_stats_ht_path
            else None,
            truth_data_ht=get_truth_ht(),
            allele_data_ht=hl.read_table(allele_data_ht_path),
            allele_counts_ht=hl.read_table(qc_ac_ht_path),
            impute_features=impute_features,
            use_adj_genotypes=use_adj_genotypes,
            n_partitions=n_partitions,
            work_bucket=work_bucket,
        )
        annotations_ht.write(out_annotations_ht_path, overwrite=True)
        logger.info(
            f'Completed annotation wrangling for random '
            f'forests model training, written to {out_annotations_ht_path}'
        )

    if not overwrite and utils.file_exists(out_results_ht_path):
        # Results already exist
        logger.info(f'Reusing {out_results_ht_path}')
    else:
        if model_file_path:
            # Overwriting, but want to reuse preexisting model
            assert model_file_path and file_exists(model_file_path)
            assert model_training_ht_path and file_exists(model_training_ht_path)
            logger.info(
                f'Loading provided model file {model_file_path}, '
                f'training ht {model_training_ht_path}'
            )
            model = load_model(model_file_path)
            training_ht = hl.read_table(model_training_ht_path)
        else:
            # Train a new model
            logger.info('Training new model')
            model, model_id, training_ht = train_model(
                annotations_ht=annotations_ht,
                rf_json_path=rf_json_path,
                model_id=model_id,
                work_bucket=work_bucket,
                overwrite=overwrite,
                use_adj_genotypes=use_adj_genotypes,
                filter_centromere_telomere=filter_centromere_telomere,
                no_transmitted_singletons=no_transmitted_singletons,
                no_inbreeding_coeff=no_inbreeding_coeff,
                fp_to_tp=fp_to_tp,
                max_depth=max_depth,
                num_trees=num_trees,
                test_intervals=test_intervals,
                vqsr_filters_split_ht_path=vqsr_filters_split_ht_path,
            )
            logger.info('Done training model')

        logger.info(f'Applying RF model...')
        ht = apply_rf_model(
            training_ht,
            rf_model=model,
            features=hl.eval(training_ht.features),
            label=LABEL_COL,
        )

        ht.write(out_results_ht_path, overwrite=True)

        summary_ht = ht.group_by(
            'tp', 'fp', TRAIN_COL, LABEL_COL, PREDICTION_COL
        ).aggregate(n=hl.agg.count())

        summary_ht.show(n=20)


def train_model(
    annotations_ht: hl.Table,
    rf_json_path: str,
    model_id: Optional[str],
    work_bucket: str,
    overwrite: bool,
    use_adj_genotypes: bool,
    filter_centromere_telomere: Optional[bool],
    no_transmitted_singletons: Optional[bool],
    no_inbreeding_coeff: Optional[bool],
    fp_to_tp: Optional[float],
    max_depth: Optional[int],
    num_trees: Optional[int],
    test_intervals: Optional[List[str]],
    vqsr_filters_split_ht_path: Optional[str],
) -> Tuple[hl.Table, str, pyspark.ml.PipelineModel]:
    """
    Train RF model
    :param annotations_ht:
    :param rf_json_path:
    :param model_id:
    :param work_bucket:
    :param overwrite:
    :param use_adj_genotypes:
    :param filter_centromere_telomere:
    :param no_transmitted_singletons:
    :param no_inbreeding_coeff:
    :param fp_to_tp:
    :param max_depth:
    :param num_trees:
    :param test_intervals:
    :param vqsr_filters_split_ht_path:
    :return: model, model_id, training_ht
    """
    rf_runs = get_rf_runs(rf_json_path)
    if not model_id:
        model_id = f'rf_{str(uuid.uuid4())[:8]}'
        while model_id in rf_runs:
            model_id = f'rf_{str(uuid.uuid4())[:8]}'
    training_ht, model = train_rf(
        annotations_ht,
        vqsr_filters_split_ht=hl.read_table(vqsr_filters_split_ht_path)
        if vqsr_filters_split_ht_path
        else None,
        fp_to_tp=fp_to_tp,
        num_trees=num_trees,
        max_depth=max_depth,
        no_transmitted_singletons=no_transmitted_singletons,
        no_inbreeding_coeff=no_inbreeding_coeff,
        filter_centromere_telomere=filter_centromere_telomere,
        test_intervals=test_intervals,
    )
    out_model_training_ht_path = join(
        work_bucket, 'models', f'{model_id}-rf_training.ht'
    )
    out_model_file_path = join(work_bucket, 'models', f'{model_id}-rf.model')

    logger.info(f'Saving training Table to {out_model_training_ht_path}')
    training_ht = training_ht.checkpoint(
        out_model_training_ht_path, overwrite=overwrite
    )
    logger.info('Adding run to RF run list')
    rf_runs[model_id] = get_run_data(
        input_args={
            'transmitted_singletons': None
            if vqsr_filters_split_ht_path
            else not no_transmitted_singletons,
            'adj': use_adj_genotypes,
            'vqsr_training': vqsr_filters_split_ht_path is not None,
            'filter_centromere_telomere': filter_centromere_telomere,
        },
        test_intervals=test_intervals if test_intervals else None,
        features_importance=dict(hl.eval(training_ht.features_importance)),
        test_results=hl.eval(training_ht.test_results),
    )
    with hl.hadoop_open(rf_json_path, 'w') as f:
        json.dump(rf_runs, f)

    logger.info(f'Saving RF model to {out_model_file_path}')
    save_model(model, out_model_file_path, overwrite=overwrite)
    return model, model_id, training_ht


def create_rf_ht(
    info_ht: hl.Table,
    inbreeding_ht: hl.Table,
    trio_stats_ht: hl.Table,
    allele_data_ht: hl.Table,
    allele_counts_ht: hl.Table,
    truth_data_ht: hl.Table,
    impute_features: bool = True,
    use_adj_genotypes: bool = False,
    n_partitions: int = 5000,
    work_bucket: str = None,
) -> hl.Table:
    """
    Creates a Table with all necessary annotations for the random forest model.
    Annotations that are included:
        Features for RF:
            - InbreedingCoeff
            - variant_type
            - allele_type
            - n_alt_alleles
            - has_star
            - AS_QD
            - AS_pab_max
            - AS_MQRankSum
            - AS_SOR
            - AS_ReadPosRankSum
        Training sites (bool):
            - transmitted_singleton
            - fail_hard_filters - (ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30)
    :param bool impute_features: Whether to impute features using feature medians (
    this is done by variant type)
    :param str use_adj_genotypes: Whether to use adj genotypes
    :param int n_partitions: Number of partitions to use for final annotated table
    :param str checkpoint_path: Optional checkpoint path for the Table before median
    imputation and/or aggregate summary
    :return: Hail Table ready for RF
    :rtype: Table
    """

    group = 'adj' if use_adj_genotypes else 'raw'

    ht = info_ht
    ht = ht.transmute(**ht.info)
    ht = ht.select('lowqual', 'AS_lowqual', 'FS', 'MQ', 'QD', *INFO_FEATURES)

    inbreeding_ht = inbreeding_ht.select(
        InbreedingCoeff=hl.if_else(
            hl.is_nan(inbreeding_ht.InbreedingCoeff),
            hl.null(hl.tfloat32),
            inbreeding_ht.InbreedingCoeff,
        )
    )
    if trio_stats_ht is not None:
        trio_stats_ht = trio_stats_ht.select(
            f'n_transmitted_{group}', f'ac_children_{group}'
        )

    logger.info('Annotating Table with all columns from multiple annotation Tables')
    ht = ht.annotate(
        **inbreeding_ht[ht.key],
        **truth_data_ht[ht.key],
        **allele_data_ht[ht.key].allele_data,
        **allele_counts_ht[ht.key],
    )
    if trio_stats_ht is not None:
        ht = ht.annotate(
            **trio_stats_ht[ht.key],
        )

    # Filter to only variants found in high quality samples and are not lowqual
    ht = ht.filter((ht[f'ac_qc_samples_{group}'] > 0) & ~ht.AS_lowqual)
    ht = ht.select(
        'a_index',
        'was_split',
        *FEATURES,
        *TRUTH_DATA,
        fail_hard_filters=(ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30),
        transmitted_singleton=ht[f'n_transmitted_{group}']
        == 1 & (ht[f'ac_qc_samples_{group}'] == 2)
        if trio_stats_ht is not None
        else hl.literal(False),
        singleton=ht.ac_release_samples_raw == 1,
        ac_raw=ht.ac_qc_samples_raw,
        ac=ht.ac_release_samples_adj,
        ac_qc_samples_unrelated_raw=ht.ac_qc_samples_unrelated_raw,
    )

    ht = ht.repartition(n_partitions, shuffle=False)
    if work_bucket:
        checkpoint_path = join(work_bucket, 'rf-annotations-before-impute.ht')
        ht = ht.checkpoint(checkpoint_path, overwrite=True)

    if impute_features:
        ht = median_impute_features(ht, {'variant_type': ht.variant_type})

    summary = ht.group_by(
        'omni',
        'mills',
        'transmitted_singleton',
    ).aggregate(n=hl.agg.count())
    logger.info('Summary of truth data annotations:')
    summary.show(20)

    return ht


def train_rf(
    ht: hl.Table,
    vqsr_filters_split_ht: Optional[hl.Table],
    fp_to_tp: Optional[float] = 1.0,
    num_trees: Optional[int] = 500,
    max_depth: Optional[int] = 5,
    no_transmitted_singletons: Optional[bool] = False,
    no_inbreeding_coeff: Optional[bool] = False,
    filter_centromere_telomere: Optional[bool] = False,
    test_intervals: Optional[List[str]] = None,
) -> Tuple[hl.Table, pyspark.ml.PipelineModel]:
    """
    Train random forest model using `train_rf_model`
    :param ht: Table containing annotations needed for RF training, built with `create_rf_ht()`
    :param vqsr_filters_split_ht: Use VQSR these training sites to train the RF. Table is build with load_vqsr.py
    :param fp_to_tp: Ratio of FPs to TPs for creating the RF model. If set to 0, all training examples are used.
    :param num_trees: Number of trees in the RF model.
    :param max_depth: Maxmimum tree depth in the RF model.
    :param no_transmitted_singletons: Do not use transmitted singletons for training.
    :param no_inbreeding_coeff: Do not use inbreeding coefficient as a feature for training.
    :param filter_centromere_telomere: Filter centromeres and telomeres before training.
    :param test_intervals: Specified interval(s) will be held out for testing and evaluation only. (default to 'chr20')
    :return: `ht` annotated with training information and the RF model
    """
    features = FEATURES
    if no_inbreeding_coeff:
        logger.info('Removing InbreedingCoeff from list of features...')
        features.remove('InbreedingCoeff')

    if vqsr_filters_split_ht is not None:
        logger.info('Using VQSR training sites for RF training...')
        ht = ht.annotate(
            vqsr_POSITIVE_TRAIN_SITE=vqsr_filters_split_ht[
                ht.key
            ].info.POSITIVE_TRAIN_SITE,
            vqsr_NEGATIVE_TRAIN_SITE=vqsr_filters_split_ht[
                ht.key
            ].info.NEGATIVE_TRAIN_SITE,
        )
        tp_expr = ht.vqsr_POSITIVE_TRAIN_SITE
        fp_expr = ht.vqsr_NEGATIVE_TRAIN_SITE
    else:
        fp_expr = ht.fail_hard_filters
        tp_expr = ht.omni | ht.mills
        if not no_transmitted_singletons:
            tp_expr = tp_expr | ht.transmitted_singleton

    if test_intervals is None:
        test_intervals = ['chr20']
    logger.info(f'Using test intervals: {test_intervals}')
    test_intervals = [
        hl.parse_locus_interval(x, reference_genome='GRCh38') for x in test_intervals
    ]

    ht = ht.annotate(tp=tp_expr, fp=fp_expr)

    if filter_centromere_telomere:
        logger.info('Filtering centromeres and telomeres from HT...')
        rf_ht = ht.filter(~hl.is_defined(telomeres_and_centromeres.ht()[ht.locus]))
    else:
        rf_ht = ht

    trained_rf_ht, rf_model = train_rf_model(
        rf_ht,
        rf_features=features,
        tp_expr=rf_ht.tp,
        fp_expr=rf_ht.fp,
        fp_to_tp=fp_to_tp,
        num_trees=num_trees,
        max_depth=max_depth,
        test_expr=hl.literal(test_intervals).any(
            lambda interval: interval.contains(rf_ht.locus)
        ),
    )
    logger.info('Joining original RF Table with training information')
    ht = ht.join(trained_rf_ht, how='left')
    return ht, rf_model


if __name__ == '__main__':
    main()  # pylint: disable=E1120
