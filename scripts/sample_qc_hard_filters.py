#!/usr/bin/env python

"""
Run sample QC on a MatrixTable, hard filter samples, add soft filter labels.
"""

import logging

import click
import hail as hl
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path
from gnomad.sample_qc.pipeline import annotate_sex
from hail.vds import read_vds, filter_chromosomes, sample_qc, filter_intervals

from larcoh import utils

logger = logging.getLogger('sample_qc')
logger.setLevel('INFO')


@click.command()
@click.option(
    '--vds',
    'vds_path',
    required=True,
    callback=utils.get_validation_callback(ext='vds', must_exist=True),
    help='path to the input dataset',
)
@click.option(
    '--cohort-tsv',
    'cohort_tsv_path',
    callback=utils.get_validation_callback(ext='tsv', must_exist=True),
    required=True,
)
@click.option(
    '--out-ht',
    'out_ht_path',
    callback=utils.get_validation_callback(ext='ht'),
    required=True,
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,missing-function-docstring
    vds_path: str,
    cohort_tsv_path: str,
    out_ht_path: str,
):
    vds = read_vds(vds_path)
    vds = filter_chromosomes(vds, keep=[f'chr{chrom}' for chrom in range(1, 23)])

    # Remove centromeres and telomeres in case they were included and any reference blocks
    tel_cent_ht = hl.read_table(str(reference_path('gnomad/tel_and_cent_ht')))
    vds = filter_intervals(vds, tel_cent_ht, keep=False)

    ht = vds.variant_data.cols()
    ht = ht.annotate(**hl.import_table(cohort_tsv_path))
    ht = ht.annotate(**sample_qc(vds))

    logger.info('Inferring sex')
    # `sex_ht` row fields: is_female, chr20_mean_dp, sex_karyotype
    ht = ht.annotate(**annotate_sex(vds, gt_expr='LGT'))

    # `out_hard_filtered_samples_ht_path` row fields: hard_filters;
    # also includes only failed samples
    ht = compute_hard_filters(ht)

    ht = ht.repartition(100)
    ht.checkpoint(out_ht_path, overwrite=True)


def compute_hard_filters(ht: hl.Table) -> hl.Table:
    """
    Uses the sex imputation results, variant sample qc, and input QC metrics
    to apply filters to the samples in vds, and create a table with
    only samples that fail at least one sample.
    """
    logger.info('Generating hard filters')

    ht = ht.annotate(hard_filters=hl.empty_set(hl.tstr))

    # Helper function to add filters into the `hard_filters` set
    def add_filter(ht_, expr, name):
        return ht_.annotate(
            hard_filters=hl.if_else(
                expr & hl.is_defined(expr), ht_.hard_filters.add(name), ht_.hard_filters
            )
        )

    # Remove samples with ambiguous sex assignments
    ht = add_filter(ht, ht.sex_karyotype == 'ambiguous', 'ambiguous_sex')
    ht = add_filter(
        ht,
        ~hl.set({'ambiguous', 'XX', 'XY'}).contains(ht.sex_karyotype),
        'sex_aneuploidy',
    )

    # Remove low-coverage samples
    # chrom 20 coverage is computed to infer sex and used here
    ht = add_filter(
        ht,
        ht.chr20_mean_dp < get_config()['larhoc']['hardfiltering']['min_coverage'],
        'low_coverage',
    )

    # Remove extreme raw bi-allelic sample QC outliers
    ht = add_filter(
        ht,
        (
            (ht.sample_qc.n_snp > get_config()['larhoc']['hardfiltering']['max_n_snps'])
            | (
                ht.sample_qc.n_snp
                < get_config()['larhoc']['hardfiltering']['min_n_snps']
            )
            | (
                ht.sample_qc.n_singleton
                > get_config()['larhoc']['hardfiltering']['max_n_singletons']
            )
            | (
                ht.sample_qc.r_het_hom_var
                > get_config()['larhoc']['hardfiltering']['max_r_het_hom']
            )
        ),
        'bad_biallelic_metrics',
    )

    ht = add_filter(
        ht,
        hl.is_defined(ht.r_contamination)
        & (
            ht.r_contamination
            > get_config()['larhoc']['hardfiltering']['max_r_contamination']
        ),
        'contamination',
    )
    ht = add_filter(
        ht,
        hl.is_defined(ht.r_duplication)
        & (
            ht.r_duplication
            > get_config()['larhoc']['hardfiltering']['max_r_duplication']
        ),
        'dup_rate',
    )
    ht = add_filter(
        ht,
        hl.is_defined(ht.insert_size)
        & (ht.insert_size < get_config()['larhoc']['hardfiltering']['min_insert_size']),
        'insert_size',
    )
    ht = ht.annotate_globals(
        hard_filter_cutoffs=hl.struct(**get_config()['larhoc']['hardfiltering'])
    )
    ht = ht.filter(hl.len(ht.hard_filters) > 0)
    return ht


if __name__ == '__main__':
    main()  # pylint: disable=E1120
