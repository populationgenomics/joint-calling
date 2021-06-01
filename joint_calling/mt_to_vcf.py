"""
Essentially a verbatim copy of: hail-ukbb-200k-callset:mt_to_vcf.py
Source:
https://github.com/populationgenomics/hail-ukbb-200k-callset/blob/master/mt_to_vcf.py
"""
import hail as hl

from gnomad.utils.vcf import ht_to_vcf_mt
from gnomad.utils.sparse_mt import default_compute_info


def mt_to_sites_only_mt(mt: hl.MatrixTable, n_partitions: int):
    """
    Convert matrix table (mt) into sites-only VCF by applying operations:
        - filter_rows_and_add_tags
        - create_info_ht
        - ht_to_vcf_mt

    :return: hl.MatrixTable
    """

    # chain these operations together
    operations = [
        filter_rows_and_add_tags,
        lambda mt: create_info_ht(mt, n_partitions=n_partitions),
        ht_to_vcf_mt,
    ]

    inter_mts = [mt]
    for op in operations:
        # perform the next operation, on the result of the previous one
        inter_mts.append(op(inter_mts[-1]))

    final_mt = inter_mts[-1]
    return final_mt


def filter_rows_and_add_tags(mt: hl.MatrixTable):
    """Filter rows and add tags"""
    mt = hl.experimental.densify(mt)
    # Filter to only non-reference sites
    mt = mt.filter_rows((hl.len(mt.alleles) > 1) & (hl.agg.any(mt.LGT.is_non_ref())))

    # annotate site level DP as site_dp onto the mt rows to avoid name collision
    mt = mt.annotate_rows(site_dp=hl.agg.sum(mt.DP))

    # Add AN tag as ANS
    return mt.annotate_rows(ANS=hl.agg.count_where(hl.is_defined(mt.LGT)) * 2)


def create_info_ht(mt: hl.MatrixTable, n_partitions: int):
    """Create info table from vcf matrix table"""
    info_ht = default_compute_info(mt, site_annotations=True, n_partitions=n_partitions)
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(DP=mt.rows()[info_ht.key].site_dp)
    )
    return info_ht
