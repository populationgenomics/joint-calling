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

from larcoh.query_utils import get_validation_callback

logger = logging.getLogger()
logger.setLevel('INFO')


@click.command()
@click.option(
    '--vds',
    'vds_path',
    required=True,
    callback=get_validation_callback(ext='vds', must_exist=True),
    help='path to the input dataset',
)
@click.option(
    '--cohort-tsv',
    'cohort_tsv_path',
    callback=get_validation_callback(ext='tsv', must_exist=True),
    required=True,
)
@click.option(
    '--out-ht',
    'out_ht_path',
    callback=get_validation_callback(ext='ht'),
    required=True,
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,missing-function-docstring
    vds_path: str,
    cohort_tsv_path: str,
    out_ht_path: str,
):
    hl.init(default_reference='GRCh38')

    vds = read_vds(vds_path)

    # Filter to autosomes:
    vds = filter_chromosomes(vds, keep=[f'chr{chrom}' for chrom in range(1, 23)])

    # Remove centromeres and telomeres:
    tel_cent_ht = hl.read_table(str(reference_path('gnomad/tel_and_cent_ht')))
    vds = filter_intervals(vds, tel_cent_ht, keep=False)

    # Initialise sample table:
    ht = hl.import_table(cohort_tsv_path)
    ht.describe()

    # Run hail sample QC stats:
    ht = ht.annotate(sample_qc=sample_qc(vds)[ht.s])
    ht.describe()

    # Load calling intervals
    calling_intervals_path = reference_path('broad/genome_calling_interval_lists')
    calling_intervals_ht = hl.import_locus_intervals(
        str(calling_intervals_path), reference_genome='GRCh38'
    )
    logger.info('Calling intervals table:')
    calling_intervals_ht.describe()

    # Infer sex (adds row fields: is_female, chr20_mean_dp, sex_karyotype)
    sex_ht = annotate_sex(
        vds,
        included_intervals=calling_intervals_ht,
        gt_expr='LGT',
    )
    logger.info('Sex table:')
    sex_ht.describe()
    ht = ht.annotate(**sex_ht[ht.s])
    ht.describe()

    # Populate ht.filters:
    ht = compute_hard_filters(ht)

    ht = ht.repartition(100)
    logger.info('Sample QC table:')
    ht.describe()
    ht.checkpoint(out_ht_path, overwrite=True)


def compute_hard_filters(ht: hl.Table) -> hl.Table:
    """
    Uses the sex imputation results, variant sample qc, and input QC metrics
    to apply filters to the samples in vds, and create a table with
    only samples that fail at least one sample.
    """
    logger.info('Generating hard filters')
    ht = ht.annotate(filters=hl.empty_set(hl.tstr))

    # Helper function to add filters into the `hard_filters` set
    def add_filter(ht_, expr, name):
        return ht_.annotate(
            filters=hl.if_else(
                expr & hl.is_defined(expr), ht_.filters.add(name), ht_.filters
            )
        )

    # Remove samples with ambiguous sex assignments
    ht = add_filter(ht, ht.sex_karyotype == 'ambiguous', 'ambiguous_sex')
    ht = add_filter(
        ht,
        ~hl.set({'ambiguous', 'XX', 'XY'}).contains(ht.sex_karyotype),
        'sex_aneuploidy',
    )

    cutoffs = get_config()['larhoc']['hardfiltering']

    # Remove low-coverage samples
    # chrom 20 coverage is computed to infer sex and used here
    ht = add_filter(
        ht,
        ht.chr20_mean_dp < cutoffs['min_coverage'],
        'low_coverage',
    )

    # Remove extreme raw bi-allelic sample QC outliers
    ht = add_filter(
        ht,
        (
            (ht.sample_qc.n_snp > cutoffs['max_n_snps'])
            | (ht.sample_qc.n_snp < cutoffs['min_n_snps'])
            | (ht.sample_qc.n_singleton > cutoffs['max_n_singletons'])
            | (ht.sample_qc.r_het_hom_var > cutoffs['max_r_het_hom'])
        ),
        'bad_biallelic_metrics',
    )

    ht = add_filter(
        ht,
        hl.is_defined(ht.r_contamination)
        & (ht.r_contamination > cutoffs['max_r_contamination']),
        'contamination',
    )
    ht = add_filter(
        ht,
        hl.is_defined(ht.r_duplication)
        & (ht.r_duplication > cutoffs['max_r_duplication']),
        'dup_rate',
    )
    ht = add_filter(
        ht,
        hl.is_defined(ht.insert_size) & (ht.insert_size < cutoffs['min_insert_size']),
        'insert_size',
    )
    ht = ht.annotate_globals(hard_filter_cutoffs=hl.struct(**cutoffs))
    return ht


if __name__ == '__main__':
    main()  # pylint: disable=E1120
