#!/usr/bin/env python
# pylint: skip-file

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

from gnomad.variant_qc.pipeline import train_rf_model
from gnomad.variant_qc.random_forest import (
    apply_rf_model,
    get_rf_runs,
    get_run_data,
    load_model,
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
    '--annotations-ht',
    'annotations_ht_path',
    required=True,
    callback=get_validation_callback(ext='ht', must_exist=True),
    help='RF annotations, generated by create_rf_annotations.py',
)
@click.option(
    '--vqsr-filters-split-ht',
    'vqsr_filters_split_ht_path',
    callback=get_validation_callback(ext='ht'),
    help='Use VQSR training sites to train the RF (generated by load_vqsr.py)',
)
@click.option(
    '--out-results-ht',
    'out_results_ht_path',
    required=True,
    callback=get_validation_callback(ext='ht'),
    help='Path to write to RF results to',
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
    '--model-file',
    'model_file_path',
    help='Pre-created model file (requires --model-id and --model-training-ht)',
)
@click.option(
    '--model-training-ht',
    'model_training_ht_path',
    help='Pre-created model training file (requires --model-id and --model-file)',
)
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
@click.option(
    '--vqsr-model-id',
    'vqsr_model_id',
    help='If a VQSR model ID is provided the VQSR training annotations will be '
    'used for training.',
    default='vqsr_alleleSpecificTrans',
    type=click.Choice(
        ['vqsr_classic', 'vqsr_alleleSpecific', 'vqsr_alleleSpecificTrans'],
        case_sensitive=False,
    ),
)
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
def main(  # pylint: disable=too-many-arguments,too-many-locals,missing-function-docstring
    annotations_ht_path: str,
    vqsr_filters_split_ht_path: Optional[str],
    out_results_ht_path: str,
    model_id: str,
    work_bucket: str,
    local_tmp_dir: str,
    overwrite: bool,
    model_file_path: Optional[str],
    model_training_ht_path: Optional[str],
    fp_to_tp: Optional[float],
    test_intervals: Optional[List[str]],
    num_trees: Optional[int],
    max_depth: Optional[int],
    use_adj_genotypes: bool,
    vqsr_model_id: Optional[str],  # pylint: disable=unused-argument
    no_transmitted_singletons: Optional[bool],
    no_inbreeding_coeff: Optional[bool],
    filter_centromere_telomere: Optional[bool],
):
    local_tmp_dir = utils.init_hail('variant_qc_random_forest', local_tmp_dir)

    rf_json_path = join(work_bucket, 'rf_runs.json')
    rf_runs = get_rf_runs(rf_json_path)
    if rf_runs:
        logger.info(f'RF runs:')
        pretty_print_runs(rf_runs)

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
            annotations_ht=hl.read_table(annotations_ht_path),
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
        rf_ht = ht.filter(
            ~hl.is_defined(hl.read_table(utils.TEL_AND_CENT_HT)[ht.locus])
        )
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
