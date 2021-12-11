"""
Utility module to encapsulate pointers to the reference files used in the pipeline
"""

import functools
import logging
import operator
from os.path import join
from typing import List, Optional, Union

import hail as hl


logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


REF_BUCKET = 'gs://cpg-reference'
NOALT_REGIONS = join(REF_BUCKET, 'noalt.bed')
SOMALIER_SITES = join(REF_BUCKET, 'somalier/v0/sites.hg38.vcf.gz')

BROAD_REF_BUCKET = f'{REF_BUCKET}/hg38/v1'
REF_FASTA = join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
UNPADDED_INTERVALS = join(BROAD_REF_BUCKET, 'hg38.even.handcurated.20k.intervals')

GNOMAD_REF_BUCKET = f'{REF_BUCKET}/gnomad/v0'
TEL_AND_CENT_HT = join(
    GNOMAD_REF_BUCKET,
    'telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht',
)
LCR_INTERVALS_HT = join(GNOMAD_REF_BUCKET, 'lcr_intervals/LCRFromHengHg38.ht')
SEG_DUP_INTERVALS_HT = join(GNOMAD_REF_BUCKET, 'seg_dup_intervals/GRCh38_segdups.ht')
CLINVAR_HT = join(GNOMAD_REF_BUCKET, 'clinvar/clinvar_20190923.ht')
HAPMAP_HT = join(GNOMAD_REF_BUCKET, 'hapmap/hapmap_3.3.hg38.ht')
KGP_OMNI_HT = join(GNOMAD_REF_BUCKET, 'kgp/1000G_omni2.5.hg38.ht')
KGP_HC_HT = join(GNOMAD_REF_BUCKET, 'kgp/1000G_phase1.snps.high_confidence.hg38.ht')
MILLS_HT = join(
    GNOMAD_REF_BUCKET, 'mills/Mills_and_1000G_gold_standard.indels.hg38.ht/'
)

GENCODE_GTF = join(REF_BUCKET, 'gencode/gencode.v29.annotation.gtf.bgz')

GNOMAD_HT = (
    'gs://gcp-public-data--gnomad/release/3.1/ht/genomes/gnomad.genomes.v3.1.sites.ht'
)

ANCESTRY_SITES = 'gs://cpg-reference/hg38/ancestry/v3/pca_sites_90k.ht'

COHORT = 'az'

"""
Subsets in gnomAD v3.1 that are broken down by their known subpops instead
of global pops in the frequency struct.
"""
COHORTS_WITH_POP_STORED_AS_SUBPOP = ['thousand-genomes', 'hgdp']

# SUBSETS can potentially contain "non_cancer", "non_neuro", etc
SUBSETS = COHORTS_WITH_POP_STORED_AS_SUBPOP


def filter_low_conf_regions(
    mt: Union[hl.MatrixTable, hl.Table],
    filter_lcr: bool = True,
    filter_segdup: bool = True,
    filter_telomeres_and_centromeres: bool = False,
    high_conf_regions: Optional[List[str]] = None,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter low-confidence regions.

    :param mt: MatrixTable or Table to filter
    :param filter_lcr: Whether to filter LCR regions
    :param filter_segdup: Whether to filter Segdup regions
    :param filter_telomeres_and_centromeres: Whether to filter telomeres and centromeres
    :param high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    :return: MatrixTable or Table with low confidence regions removed
    """
    criteria = []
    if filter_lcr:
        lcr = hl.read_table(LCR_INTERVALS_HT)
        criteria.append(hl.is_missing(lcr[mt.locus]))

    if filter_segdup:
        segdup = hl.read_table(SEG_DUP_INTERVALS_HT)
        criteria.append(hl.is_missing(segdup[mt.locus]))

    if filter_telomeres_and_centromeres:
        telomeres_and_centromeres = hl.read_table(TEL_AND_CENT_HT)
        criteria.append(hl.is_missing(telomeres_and_centromeres[mt.locus]))

    if high_conf_regions is not None:
        for region in high_conf_regions:
            region = hl.import_locus_intervals(region)
            criteria.append(hl.is_defined(region[mt.locus]))

    if criteria:
        filter_criteria = functools.reduce(operator.iand, criteria)
        if isinstance(mt, hl.MatrixTable):
            mt = mt.filter_rows(filter_criteria)
        else:
            mt = mt.filter(filter_criteria)

    return mt


def get_truth_ht() -> hl.Table:
    """
    Return a table with annotations from the latest version of the corresponding truth data.

    The following annotations are included:
        - hapmap
        - kgp_omni (1000 Genomes intersection Onni 2.5M array)
        - kgp_phase_1_hc (high confidence sites in 1000 genonmes)
        - mills (Mills & Devine indels)

    :return: A table with the latest version of popular truth data annotations
    """
    return (
        hl.read_table(HAPMAP_HT)
        .select(hapmap=True)
        .join(hl.read_table(KGP_OMNI_HT).select(omni=True), how='outer')
        .join(hl.read_table(KGP_HC_HT).select(kgp_phase1_hc=True), how='outer')
        .join(hl.read_table(MILLS_HT).select(mills=True), how='outer')
        .repartition(200, shuffle=False)
        .persist()
    )
