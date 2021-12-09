#!/usr/bin/env python

"""
Computes a HT with the typical GATK AS and site-level info fields
as well as ACs and lowqual fields. Note that this table doesn't 
split multi-allelic sites.

Generates info.ht, info-split.ht, needed for sample_qc and for 
random forest jobs.
"""

import logging
import click
import hail as hl

from gnomad.utils.annotations import (
    get_adj_expr,
    get_lowqual_expr,
)
from gnomad.utils.sparse_mt import (
    get_as_info_expr,
    get_site_info_expr,
    INFO_INT32_SUM_AGG_FIELDS,
    INFO_SUM_AGG_FIELDS,
    split_info_annotation,
    split_lowqual_annotation,
)

from joint_calling.utils import file_exists
from joint_calling import utils, _version

logger = logging.getLogger('qc-annotations')
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--out-info-ht',
    'out_info_ht_path',
    required=True,
)
@click.option(
    '--out-split-info-ht',
    'out_split_info_ht_path',
    required=True,
)
@click.option(
    '--mt',
    'mt_path',
    required=True,
    callback=utils.get_validation_callback(ext='mt', must_exist=True),
    help='path to the input MatrixTable',
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
def main(  # pylint: disable=too-many-arguments,too-many-locals,too-many-statements,missing-function-docstring
    out_info_ht_path: str,
    out_split_info_ht_path: str,
    mt_path: str,
    local_tmp_dir: str,
    overwrite: bool,
):
    local_tmp_dir = utils.init_hail('generate_info_ht', local_tmp_dir)

    all_samples_mt = utils.get_mt(mt_path)

    compute_info(
        mt=all_samples_mt,
        out_ht_path=out_info_ht_path,
        out_split_ht_path=out_split_info_ht_path,
        overwrite=overwrite,
    )


def compute_info(
    mt: hl.MatrixTable,
    out_ht_path: str,
    out_split_ht_path: str,
    overwrite: bool = False,
) -> hl.Table:
    """
    Computes a HT with the typical GATK AS and site-level info fields as well as ACs
    and lowqual fields.

    :param mt: full matrix table
    :param out_ht_path: where to write the info Table
    :param out_split_ht_path: if provided, in the info Table multiallelics will be split
    and the Table will be written to this file
    :param overwrite: overwrite checkpoints if they exist
    :return: Table with info fields
    """
    if not overwrite and file_exists(out_ht_path):
        info_ht = hl.read_table(out_ht_path)
    else:
        mt = mt.filter_rows((hl.len(mt.alleles) > 1))
        mt = mt.transmute_entries(**mt.gvcf_info)
        mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))

        # Compute AS and site level info expr
        # Note that production defaults have changed:
        # For new releases, the `RAW_MQandDP` field replaces the `RAW_MQ` and `MQ_DP` fields
        info_expr = get_site_info_expr(
            mt,
            sum_agg_fields=INFO_SUM_AGG_FIELDS,
            int32_sum_agg_fields=INFO_INT32_SUM_AGG_FIELDS,
            array_sum_agg_fields=['SB', 'RAW_MQandDP'],
        )
        info_expr = info_expr.annotate(
            **get_as_info_expr(
                mt,
                sum_agg_fields=INFO_SUM_AGG_FIELDS,
                int32_sum_agg_fields=INFO_INT32_SUM_AGG_FIELDS,
                array_sum_agg_fields=['SB', 'RAW_MQandDP'],
            )
        )

        # Add AC and AC_raw:
        # First compute ACs for each non-ref allele, grouped by adj
        grp_ac_expr = hl.agg.array_agg(
            lambda ai: hl.agg.filter(
                mt.LA.contains(ai),
                hl.agg.group_by(
                    get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD),
                    hl.agg.sum(
                        mt.LGT.one_hot_alleles(mt.LA.map(hl.str))[mt.LA.index(ai)]
                    ),
                ),
            ),
            mt.alt_alleles_range_array,
        )

        # Then, for each non-ref allele, compute
        # AC as the adj group
        # AC_raw as the sum of adj and non-adj groups
        info_expr = info_expr.annotate(
            AC_raw=grp_ac_expr.map(
                lambda i: hl.int32(i.get(True, 0) + i.get(False, 0))
            ),
            AC=grp_ac_expr.map(lambda i: hl.int32(i.get(True, 0))),
        )

        # Annotating raw MT with pab max
        info_expr = info_expr.annotate(
            AS_pab_max=hl.agg.array_agg(
                lambda ai: hl.agg.filter(
                    mt.LA.contains(ai) & mt.LGT.is_het(),
                    hl.agg.max(
                        hl.binom_test(
                            mt.LAD[mt.LA.index(ai)], hl.sum(mt.LAD), 0.5, 'two-sided'
                        )
                    ),
                ),
                mt.alt_alleles_range_array,
            )
        )

        info_ht = mt.select_rows(info=info_expr).rows()

        # Add lowqual flag
        info_ht = info_ht.annotate(
            lowqual=get_lowqual_expr(
                info_ht.alleles,
                info_ht.info.QUALapprox,
                # The indel het prior used for gnomad v3 was 1/10k bases (phred=40).
                # This value is usually 1/8k bases (phred=39).
                indel_phred_het_prior=40,
            ),
            AS_lowqual=get_lowqual_expr(
                info_ht.alleles, info_ht.info.AS_QUALapprox, indel_phred_het_prior=40
            ),
        )
        info_ht = info_ht.naive_coalesce(7500)
        info_ht.write(out_ht_path, overwrite=True)

    if out_split_ht_path and (overwrite or not file_exists(out_split_ht_path)):
        split_info_ht = split_multiallelic_in_info_table(info_ht)
        split_info_ht.write(out_split_ht_path, overwrite=True)

    return info_ht


def split_multiallelic_in_info_table(info_ht: hl.Table) -> hl.Table:
    """
    Generates an info table that splits multi-allelic sites from the multi-allelic
    info table.
    :return: Info table with split multi-allelics
    :rtype: Table
    """

    # Create split version
    info_ht = hl.split_multi(info_ht)

    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            **split_info_annotation(info_ht.info, info_ht.a_index),
        ),
        AS_lowqual=split_lowqual_annotation(info_ht.AS_lowqual, info_ht.a_index),
    )
    return info_ht


if __name__ == '__main__':
    main()  # pylint: disable=E1120
