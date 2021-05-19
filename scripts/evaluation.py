#!/usr/bin/env python

"""
Run Random Forest variant QC
"""
from os.path import join
from typing import Optional
import logging
from pprint import pformat
import click
import hail as hl

from gnomad.resources.grch38.reference_data import clinvar, telomeres_and_centromeres
from gnomad.utils.filtering import filter_low_conf_regions, filter_to_clinvar_pathogenic
from gnomad.variant_qc.evaluation import (
    compute_binned_truth_sample_concordance,
    compute_grouped_binned_ht,
    create_truth_sample_ht,
)
from gnomad.utils.annotations import (
    annotate_adj,
)
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.sample_qc.relatedness import generate_trio_stats_expr
from gnomad.variant_qc.pipeline import create_binned_ht, score_bin_agg
from gnomad.resources.grch38 import (
    na12878_giab,
    na12878_giab_hc_intervals,
    syndip,
    syndip_hc_intervals,
)

from joint_calling.utils import get_validation_callback
from joint_calling import utils
from joint_calling import _version

logger = logging.getLogger('random_forest')
logger.setLevel('INFO')


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--info-split-ht',
    'info_split_ht_path',
    required=True,
    callback=get_validation_callback(ext='ht', must_exist=True),
    help='path to info Table with split multiallelics '
    '(generated by generate_info_ht.py --out-split-info-ht)',
)
@click.option(
    '--rf-results-ht',
    'rf_results_ht_path',
    callback=get_validation_callback(ext='ht'),
    help='RF results table (created by random_forest.py)',
)
@click.option(
    '--rf-annotations-ht',
    'rf_annotations_ht_path',
    callback=get_validation_callback(ext='ht'),
    help='RF annotations table (created by random_forest.py)',
)
@click.option(
    '--mt',
    'mt_path',
    callback=get_validation_callback(ext='mt'),
    help='Path to the matrix table',
)
@click.option(
    '--fam-stats-ht',
    'fam_stats_ht_path',
    callback=get_validation_callback(ext='ht'),
    help='path to a Table with trio stats (generated by generate_rf_annotations.py)',
)
@click.option(
    '--vqsr-filters-split-ht',
    'vqsr_filters_split_ht_path',
    callback=get_validation_callback(ext='ht'),
    help='Use VQSR training sites to train the RF (generated by load_vqsr.py)',
)
@click.option(
    '--out-bin-ht',
    'out_bin_ht',
    required=True,
    help='When set, creates file annotated with bin based on rank of VQSR/RF score.',
)
@click.option(
    '--out-aggregated-bin-ht',
    'out_aggregated_bin_ht',
    help='When set, creates a file with aggregate counts of variants based on bins.',
)
@click.option(
    '--run-sanity-checks',
    'run_sanity_checks',
    is_flag=True,
    help='When set, runs ranking sanity checks.',
)
@click.option(
    '--n-bins',
    'n_bins',
    help='Number of bins for the binned file (default: 100).',
    default=100,
    type=click.INT,
)
@click.option(
    '--n-partitions',
    'n_partitions',
    type=click.INT,
    help='Desired number of partitions for output Table/MatrixTable.',
    default=5000,
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
def main(  # pylint: disable=too-many-arguments,too-many-locals
    info_split_ht_path: str,
    mt_path: str,
    fam_stats_ht_path: str,
    rf_results_ht_path: Optional[str],
    vqsr_filters_split_ht_path: Optional[str],
    rf_annotations_ht_path: Optional[str],
    out_bin_ht: str,
    out_aggregated_bin_ht: str,
    run_sanity_checks: bool,
    n_bins: int,
    n_partitions: int,
    work_bucket: str,
    local_tmp_dir: str,
    overwrite: bool,
):
    """
    Evaluate random forest and VQSR results against the truth data
    """
    local_tmp_dir = utils.init_hail('variant_qc_evaluate', local_tmp_dir)

    if overwrite or not utils.file_exists(out_bin_ht):
        scores_ht = create_bin_ht(
            info_split_ht=hl.read_table(info_split_ht_path),
            n_bins=n_bins,
            rf_results_ht=hl.read_table(rf_results_ht_path)
            if rf_results_ht_path
            else None,
            vqsr_filters_split_ht=hl.read_table(vqsr_filters_split_ht_path)
            if vqsr_filters_split_ht_path
            else None,
            rf_annotations_ht=hl.read_table(rf_annotations_ht_path)
            if rf_annotations_ht_path
            else None,
        )
        scores_ht.write(out_bin_ht, overwrite=True)
    else:
        scores_ht = hl.read_table(out_bin_ht)

    if run_sanity_checks:
        logger.info('Running sanity checks...')
        ht = scores_ht
        logger.info(
            ht.aggregate(
                hl.struct(
                    was_biallelic=hl.agg.counter(~ht.was_split),
                    has_biallelic_rank=hl.agg.counter(hl.is_defined(ht.biallelic_bin)),
                    was_singleton=hl.agg.counter(ht.singleton),
                    has_singleton_rank=hl.agg.counter(hl.is_defined(ht.singleton_bin)),
                    was_biallelic_singleton=hl.agg.counter(
                        ht.singleton & ~ht.was_split
                    ),
                    has_biallelic_singleton_rank=hl.agg.counter(
                        hl.is_defined(ht.biallelic_singleton_bin)
                    ),
                )
            )
        )

    if not fam_stats_ht_path or not utils.file_exists(fam_stats_ht_path):
        hard_filtered_samples_ht_path = (
            'gs://cpg-tob-wgs-temporary/joint_vcf/v1/work/combiner/hard_filters.ht/'
        )
        meta_ht_path = 'gs://cpg-tob-wgs-temporary/joint_vcf/v1/work/combiner/meta.ht/'
        hard_filtered_mt = utils.get_mt(
            mt_path,
            hard_filtered_samples_to_remove_ht=hl.read_table(
                hard_filtered_samples_ht_path
            ),
            meta_ht=hl.read_table(meta_ht_path),
            add_meta=True,
        )

        trios_fam_ped_file = _make_fam_file(
            sex_ht=hl.read_table(meta_ht_path),
            work_bucket='gs://cpg-tob-wgs-temporary/joint_vcf/v1/work/variant_qc/',
        )

        generate_fam_stats(
            mt=hard_filtered_mt,
            out_fam_stats_ht_path=fam_stats_ht_path,
            trios_fam_ped_file=trios_fam_ped_file,
            overwrite=False,
        )

    if out_aggregated_bin_ht:
        if overwrite or not utils.file_exists(out_aggregated_bin_ht):
            logger.warning('Use only workers, it typically crashes with preemptibles')
            agg_ht = create_aggregated_bin_ht(
                ht=scores_ht,
                trio_stats_ht=hl.read_table(fam_stats_ht_path),
                work_bucket=work_bucket,
            )
            agg_ht.write(out_aggregated_bin_ht, overwrite=True)

    mt = utils.get_mt(mt_path)
    if all(truth_sample in mt.s.collect() for truth_sample in utils.TRUTH_GVCFS):
        _truth_concordance(
            mt,
            overwrite,
            work_bucket,
            n_partitions,
            scores_ht,
            info_split_ht_path,
            n_bins,
        )


def _truth_concordance(
    mt,
    overwrite,
    work_bucket,
    n_partitions,
    scores_ht,
    info_split_ht_path,
    n_bins,
):
    logger.info(f'Extracting truth samples from MT...')
    truth_dict = {
        utils.TRUTH_GVCFS['syndip']['name']: {
            's': utils.TRUTH_GVCFS['syndip']['name'],
            'truth_mt': syndip.mt(),
            'hc_intervals': syndip_hc_intervals.ht(),
            'mt': None,
            'ht': None,
        },
        utils.TRUTH_GVCFS['na12878']['name']: {
            's': utils.TRUTH_GVCFS['na12878']['name'],
            'truth_mt': na12878_giab.mt(),
            'hc_intervals': na12878_giab_hc_intervals.ht(),
            'mt': None,
            'ht': None,
        },
    }

    mt = mt.filter_cols(
        hl.literal([v['s'] for k, v in truth_dict.items()]).contains(mt.s)
    )
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    # Checkpoint to prevent needing to go through the large table a second time
    checkpoint_mt_path = join(work_bucket, 'tmp', 'genomes_split.mt')
    logger.info(f'Saving checkpoint to {checkpoint_mt_path}')
    mt = mt.checkpoint(checkpoint_mt_path, overwrite=overwrite)

    for truth_sample in truth_dict:
        truth_samples_mt_path = join(work_bucket, 'truth_samples', f'{truth_sample}.mt')
        if not overwrite and utils.file_exists(truth_samples_mt_path):
            truth_dict[truth_sample]['mt'] = hl.read_matrix_table(truth_samples_mt_path)
        else:
            called_truth_mt = mt.filter_cols(mt.s == truth_dict[truth_sample]['s'])
            # Filter to variants in truth data
            called_truth_mt = called_truth_mt.filter_rows(
                hl.agg.any(called_truth_mt.GT.is_non_ref())
            )
            logger.info(
                f'Saving {truth_sample} called truth sample data to '
                f'{truth_samples_mt_path}'
            )
            called_truth_mt = called_truth_mt.naive_coalesce(n_partitions)
            called_truth_mt.write(truth_samples_mt_path, overwrite=True)
            truth_dict[truth_sample]['mt'] = called_truth_mt

    for truth_sample in truth_dict:
        # Merging with truth data. Computes a table for each truth sample comparing
        # the truth sample in the callset vs the truth.
        truth_ht_path = join(work_bucket, 'truth_samples', f'{truth_sample}.ht')
        if not overwrite and utils.file_exists(truth_ht_path):
            truth_dict[truth_sample]['ht'] = hl.read_table(truth_ht_path)
        else:
            logger.info(
                f'Creating a merged table with callset truth sample and truth data '
                f'for {truth_sample}...'
            )

            # Load truth data
            mt = truth_dict[truth_sample]['mt']
            truth_hc_intervals = truth_dict[truth_sample]['hc_intervals']
            truth_mt = truth_dict[truth_sample]['truth_mt']
            truth_mt = truth_mt.key_cols_by(s=hl.str(truth_dict[truth_sample]['s']))

            # Remove low quality sites
            info_ht = hl.read_table(info_split_ht_path)
            mt = mt.filter_rows(~info_ht[mt.row_key].AS_lowqual)

            ht = create_truth_sample_ht(mt, truth_mt, truth_hc_intervals)
            ht.write(truth_ht_path, overwrite=True)
            truth_dict[truth_sample]['ht'] = ht

        # Bin truth sample concordance. Merges concordance results (callset vs.
        # truth) for a given truth sample with bins from specified model
        logger.info(f'Creating binned concordance table for {truth_sample}')
        info_ht = hl.read_table(info_split_ht_path)
        ht = truth_dict[truth_sample]['ht']
        ht = ht.filter(
            ~info_ht[ht.key].AS_lowqual
            & ~hl.is_defined(telomeres_and_centromeres.ht()[ht.locus])
        )

        logger.info('Filtering out low confidence regions and segdups...')
        ht = filter_low_conf_regions(
            ht,
            filter_lcr=True,
            filter_decoy=False,  # Set if having decoy path
            filter_segdup=True,
        )

        logger.info(
            'Loading HT containing RF or VQSR scores annotated with a bin based '
            'on the rank of score...'
        )
        metric_ht = scores_ht
        ht = ht.filter(hl.is_defined(metric_ht[ht.key]))

        ht = ht.annotate(score=metric_ht[ht.key].score)

        ht = compute_binned_truth_sample_concordance(ht, metric_ht, n_bins)
        binned_concordance_ht_path = join(
            work_bucket,
            'binned_concordance',
            f'{truth_sample}_binned_concordance.ht',
        )
        ht.write(binned_concordance_ht_path, overwrite=True)


def create_bin_ht(
    info_split_ht: hl.Table,
    n_bins: int,
    rf_results_ht: Optional[hl.Table] = None,
    vqsr_filters_split_ht: Optional[hl.Table] = None,
    rf_annotations_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Creates a table with bin annotations added for a RF or VQSR run and writes it
    to its correct location in annotations.

    :param info_split_ht: table generated by generate_info_ht.py
    :param n_bins: Number of bins to bin the data into
    :param rf_results_ht: table generated by random_forest.py
    :param vqsr_filters_split_ht: table generated by load_vqsr.py
    :param rf_annotations_ht: table generated by random_forest.py
    :return: Table with bin annotations
    """
    logger.info(f'Annotating model HT with bins using {n_bins} bins')
    if vqsr_filters_split_ht:
        logger.info(f'Using a VQSR model')

        ht = vqsr_filters_split_ht

        if rf_annotations_ht is not None:
            ht = ht.annotate(**rf_annotations_ht[ht.key])

        ht = ht.annotate(
            info=info_split_ht[ht.key].info,
            score=ht.info.AS_VQSLOD,
            positive_train_site=ht.info.POSITIVE_TRAIN_SITE,
            negative_train_site=ht.info.NEGATIVE_TRAIN_SITE,
            AS_culprit=ht.info.AS_culprit,
        )

        # Remove all samples with an undefined ac_raw because ac_raw was calculated
        # on the high quality samples only and VQSR was run before sample filtering
        ht = ht.filter(hl.is_defined(ht.ac_raw))

    else:
        logger.info(f'Using an RF model')
        ht = rf_results_ht
        ht = ht.annotate(
            info=info_split_ht[ht.key].info,
            positive_train_site=ht.tp,
            negative_train_site=ht.fp,
            score=ht.rf_probability['TP'],
        )

    ht = ht.filter(
        ~info_split_ht[ht.key].AS_lowqual
        & ~hl.is_defined(telomeres_and_centromeres.ht()[ht.locus])
    )
    ht_non_lcr = filter_low_conf_regions(
        ht,
        filter_lcr=True,
        filter_decoy=False,  # Set if having decoy path
        filter_segdup=True,
    )
    ht = ht.annotate(non_lcr=hl.is_defined(ht_non_lcr[ht.key]))
    bin_ht = create_binned_ht(ht, n_bins, add_substrat={'non_lcr': ht.non_lcr})

    return bin_ht


def create_aggregated_bin_ht(
    ht: hl.Table,
    trio_stats_ht: hl.Table,
    work_bucket: str,
) -> hl.Table:
    """
    Aggregates variants into bins, grouped by `bin_id` (rank, bi-allelic, etc.),
    contig, and `snv`, `bi_allelic`, and `singleton` status, using previously
    annotated bin information.

    For each bin, aggregates statistics needed for evaluation plots.
    :param ht: table with score bins
    :param trio_stats_ht: HT generated from a FAM file
    :param work_bucket: bucket to write temporary files to
    :return: Table of aggregate statistics by bin
    """

    # Count variants for ranking
    count_expr = {
        x: hl.agg.filter(
            hl.is_defined(ht[x]),
            hl.agg.counter(
                hl.cond(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv', 'indel')
            ),
        )
        for x in ht.row
        if x.endswith('bin')
    }
    bin_variant_counts = ht.aggregate(hl.struct(**count_expr))
    logger.info(f'Found the following variant counts:\n {pformat(bin_variant_counts)}')
    ht = ht.annotate_globals(bin_variant_counts=bin_variant_counts)

    # Load ClinVar pathogenic data
    clinvar_pathogenic_ht = filter_to_clinvar_pathogenic(clinvar.ht())
    ht = ht.annotate(clinvar_path=hl.is_defined(clinvar_pathogenic_ht[ht.key]))

    logger.info(f'Creating grouped bin table...')
    checkpoint_path = join(work_bucket, 'tmp', f'grouped_bin.ht')
    grouped_binned_ht = compute_grouped_binned_ht(ht, checkpoint_path=checkpoint_path)

    logger.info(f'Aggregating grouped bin table...')
    parent_ht = grouped_binned_ht._parent  # pylint: disable=protected-access
    agg_ht = grouped_binned_ht.aggregate(
        n_clinvar_path=hl.agg.count_where(parent_ht.clinvar_path),
        **score_bin_agg(grouped_binned_ht, fam_stats_ht=trio_stats_ht),
    )
    return agg_ht


def generate_fam_stats(
    mt: hl.MatrixTable,
    out_fam_stats_ht_path: str,
    overwrite: bool,
    trios_fam_ped_file: Optional[str] = None,
) -> hl.Table:
    """
    Calculate transmission and de novo mutation statistics using trios in the dataset.
    :param mt: input MatrixTable
    :param out_fam_stats_ht_path: path to write the resulting Table to
    :param overwrite: overwrite existing intermediate and output files
    :param trios_fam_ped_file: path to text file containing trio pedigree
    :return: Table containing trio stats
    """
    logger.info('Generate FAM stats')
    if not overwrite and utils.file_exists(out_fam_stats_ht_path):
        return hl.read_table(out_fam_stats_ht_path)

    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    # Load Pedigree data and filter MT to samples present in any of the trios
    ped = hl.Pedigree.read(trios_fam_ped_file, delimiter='\t')
    fam_ht = hl.import_fam(trios_fam_ped_file, delimiter='\t')
    fam_ht = fam_ht.annotate(fam_members=[fam_ht.id, fam_ht.pat_id, fam_ht.mat_id])
    fam_ht = fam_ht.explode('fam_members', name='s')
    fam_ht = fam_ht.key_by('s').select().distinct()

    mt = mt.filter_cols(hl.is_defined(fam_ht[mt.col_key]))
    logger.info(
        f'Generating family stats using {mt.count_cols()} samples from {len(ped.trios)} '
        f'trios.'
    )

    mt = filter_to_autosomes(mt)
    mt = annotate_adj(mt)
    mt = mt.select_entries('GT', 'GQ', 'AD', 'END', 'adj')
    mt = hl.experimental.densify(mt)
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mt = hl.trio_matrix(mt, pedigree=ped, complete_trios=True)
    trio_adj = mt.proband_entry.adj & mt.father_entry.adj & mt.mother_entry.adj

    ht = mt.select_rows(
        **generate_trio_stats_expr(
            mt,
            transmitted_strata={'raw': True, 'adj': trio_adj},
            de_novo_strata={
                'raw': True,
                'adj': trio_adj,
            },
            proband_is_female_expr=mt.is_female,
        )
    ).rows()

    ht = ht.filter(
        ht.n_de_novos_raw + ht.n_transmitted_raw + ht.n_untransmitted_raw > 0
    )
    ht = ht.repartition(10000, shuffle=False)
    ht.write(out_fam_stats_ht_path, overwrite=True)
    return ht


def _make_fam_file(sex_ht: hl.Table, work_bucket: str) -> str:
    trios_fam_ped_file = join(work_bucket, 'pedigree.fam')
    ped = hl.Pedigree(
        trios=[
            hl.Trio(s, is_female=is_female)
            for s, is_female in zip(sex_ht.s.collect(), sex_ht.is_female.collect())
        ]
    )
    ped.write(trios_fam_ped_file)
    return trios_fam_ped_file


if __name__ == '__main__':
    main()  # pylint: disable=E1120