import logging
import pickle
from typing import Dict, List, Optional, Set, Union

import click
import hail as hl

from gnomad.resources.grch38.gnomad import (
    COHORTS_WITH_POP_STORED_AS_SUBPOP,
    HGDP_POPS,
    TGP_POPS,
    TGP_POP_NAMES,
    POPS,
    SEXES,
    SUBSETS,
)
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.utils.filtering import remove_fields_from_constant
from gnomad.utils.vcf import (
    add_as_info_dict,
    adjust_vcf_incompatible_types,
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    AS_VQSR_FIELDS,
    create_label_groups,
    ENTRIES,
    FAF_POPS,
    FORMAT_DICT,
    HISTS,
    INFO_DICT,
    IN_SILICO_ANNOTATIONS_INFO_DICT,
    make_info_dict,
    make_vcf_filter_dict,
    REGION_FLAG_FIELDS,
    RF_FIELDS,
    SITE_FIELDS,
    VQSR_FIELDS,
)
from gnomad.variant_qc.pipeline import INBREEDING_COEFF_HARD_CUTOFF
from gnomad.utils.annotations import region_flag_expr

from joint_calling import utils, _version, resources
from joint_calling.utils import get_validation_callback


logger = logging.getLogger(__file__)
logger.setLevel('INFO')


# Add new site fields
# NEW_SITE_FIELDS = [
#     'monoallelic',
#     'transmitted_singleton',
# ]
# SITE_FIELDS.extend(NEW_SITE_FIELDS)

# Remove original alleles for containing non-releasable alleles
MISSING_ALLELE_TYPE_FIELDS = ['original_alleles', 'has_star']
ALLELE_TYPE_FIELDS = remove_fields_from_constant(
    ALLELE_TYPE_FIELDS, MISSING_ALLELE_TYPE_FIELDS
)

# Remove SB (not included in VCF) and SOR (doesn't exist in v3.1) from site fields
MISSING_SITES_FIELDS = ['SOR', 'SB']
SITE_FIELDS = remove_fields_from_constant(SITE_FIELDS, MISSING_SITES_FIELDS)

# Remove AS_VarDP from AS fields
MISSING_AS_FIELDS = ['AS_VarDP']
AS_FIELDS = remove_fields_from_constant(AS_FIELDS, MISSING_AS_FIELDS)

# Make subset list (used in properly filling out VCF header descriptions and naming VCF info fields)
SUBSET_LIST_FOR_VCF = SUBSETS.copy()
SUBSET_LIST_FOR_VCF.append('')

# Remove cohorts that have subpop frequencies stored as pop frequencies
# Inclusion of these subsets significantly increases the size of storage in the VCFs because of the many subpops
SUBSET_LIST_FOR_VCF = remove_fields_from_constant(
    SUBSET_LIST_FOR_VCF, COHORTS_WITH_POP_STORED_AS_SUBPOP
)

# Remove decoy from region field flag
MISSING_REGION_FIELDS = ['decoy']
REGION_FLAG_FIELDS = remove_fields_from_constant(
    REGION_FLAG_FIELDS, MISSING_REGION_FIELDS
)

# All missing fields to remove from vcf info dict
MISSING_INFO_FIELDS = (
    MISSING_ALLELE_TYPE_FIELDS
    + MISSING_AS_FIELDS
    + MISSING_REGION_FIELDS
    + MISSING_SITES_FIELDS
    + RF_FIELDS
)

# Remove unnecessary pop names from POP_NAMES dict
POPS = {pop: POP_NAMES[pop] for pop in POPS}

# Remove unnecessary pop names from FAF_POPS dict
FAF_POPS = {pop: POP_NAMES[pop] for pop in FAF_POPS}

# Get HGDP + TGP(KG) subset pop names
HGDP_TGP_KEEP_POPS = TGP_POPS + HGDP_POPS
HGDP_TGP_POPS = {}
for pop in HGDP_TGP_KEEP_POPS:
    if pop in TGP_POP_NAMES:
        HGDP_TGP_POPS[pop] = TGP_POP_NAMES[pop]
    else:
        HGDP_TGP_POPS[pop] = pop.capitalize()

# Used for HGDP + TGP subset MT VCF output only
FORMAT_DICT.update(
    {
        'RGQ': {
            'Number': '1',
            'Type': 'Integer',
            'Description': 'Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)',
        }
    }
)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--mt',
    'mt_path',
    required=True,
    callback=get_validation_callback(ext='mt', must_exist=True),
    help='path to the raw sparse Matrix Table generated by combine_gvcfs.py',
)
@click.option(
    '--out-ht',
    'out_ht_path',
    required=True,
    callback=get_validation_callback(ext='ht'),
    help='path to write Hail Table',
)
@click.option(
    '--out-vcf-header-txt',
    'out_vcf_header_txt_path',
    required=True,
    callback=get_validation_callback(ext='txt'),
)
@click.option(
    '--public-subset',
    'is_public_subset',
    is_flag=True,
    help='create a subset',
)
@click.option(
    '--test',
    'is_test',
    is_flag=True,
    required=True,
    help='Create release files using only 2 partitions on chr20, chrX, '
         'and chrY for testing purposes'
)
@click.option(
    '--local-tmp-dir',
    'local_tmp_dir',
    help='local directory for temporary files and Hail logs (must be local).',
)
def main(
    mt_path: str,
    out_ht_path: str,
    out_vcf_header_txt_path: str,
    is_test: bool,
    is_public_subset: bool,
    local_tmp_dir: str,
):  # pylint: disable=missing-function-docstring
    utils.init_hail(__file__, local_tmp_dir)
    
    mt = hl.read_matrix_table(mt_path)
    ht = mt.rows()
    if is_test:
        ht = filter_to_test(ht)

    # Setup of parameters and Table/MatrixTable
    parameter_dict = _build_parameter_dict(ht, is_public_subset)

    vcf_ht = _prepare_vcf_ht(ht, is_subset=is_public_subset)
    
    if is_public_subset:
        logger.info(
            'Loading subset MT and annotating with the prepared VCF HT for VCF export...'
        )
        entries_ht = mt.select_rows().select_entries(*ENTRIES)
        vcf_ht = ht.annotate_rows(**vcf_ht[entries_ht.row_key])

    logger.info('Cleaning up the VCF HT for final export...')
    vcf_ht = _cleanup_ht_for_vcf_export(vcf_ht)

    vcf_ht = vcf_ht.checkpoint(out_ht_path, overwrite=True)
    vcf_ht.describe()

    _prepare_vcf_header_dict(
        ht=ht,
        vcf_ht=vcf_ht,
        vcf_header_txt_path=out_vcf_header_txt_path,
        parameter_dict=parameter_dict,
        is_public_subset=is_public_subset,
    )


def _prepare_vcf_header_dict(
    ht,
    vcf_ht,
    vcf_header_txt_path,
    parameter_dict,
    is_public_subset,
):
    logger.info('Making histogram bin edges...')
    header_dict = prepare_vcf_header_dict(
        vcf_ht,
        subset_list=parameter_dict['subsets'],
        pops=parameter_dict['pops'],
        filtering_model_field=parameter_dict['filtering_model_field'],
        inbreeding_coeff_cutoff=ht.inbreeding_coeff_cutoff,
    )
    if not is_public_subset:
        header_dict.pop('format')

    logger.info('Saving header dict to pickle...')
    with hl.hadoop_open(vcf_header_txt_path, 'wb') as p:
        pickle.dump(header_dict, p, protocol=pickle.HIGHEST_PROTOCOL)


def populate_subset_info_dict(
    subset: str,
    description_text: str,
    pops: Dict[str, str] = POPS,
    faf_pops: Dict[str, str] = FAF_POPS,
    sexes: List[str] = SEXES,
    label_delimiter: str = '_',
) -> Dict[str, Dict[str, str]]:
    """
    Call `make_info_dict` to populate INFO dictionary with specific sexes, population names, and filtering allele
    frequency (faf) pops for the requested subset.

    :param subset: Sample subset in dataset.
    :param description_text: Text describing the sample subset that should be added to the INFO description.
    :param pops: Dict of sample global population names for gnomAD genomes. Default is POPS.
    :param faf_pops: Dict with faf pop names (keys) and descriptions (values).  Default is FAF_POPS.
    :param sexes: gnomAD sample sexes used in VCF export. Default is SEXES.
    :param label_delimiter: String to use as delimiter when making group label combinations. Default is '_'.
    :return: Dictionary containing Subset specific INFO header fields.
    """

    vcf_info_dict = {}
    faf_label_groups = create_label_groups(pops=faf_pops, sexes=sexes)
    for label_group in faf_label_groups:
        vcf_info_dict.update(
            make_info_dict(
                prefix=subset,
                prefix_before_metric=True if 'gnomad' in subset else False,
                pop_names=faf_pops,
                label_groups=label_group,
                label_delimiter=label_delimiter,
                faf=True,
                description_text=description_text,
            )
        )

    label_groups = create_label_groups(pops=pops, sexes=sexes)
    for label_group in label_groups:
        vcf_info_dict.update(
            make_info_dict(
                prefix=subset,
                prefix_before_metric=True if 'gnomad' in subset else False,
                pop_names=pops,
                label_groups=label_group,
                label_delimiter=label_delimiter,
                description_text=description_text,
            )
        )

    # Add popmax to info dict
    vcf_info_dict.update(
        make_info_dict(
            prefix=subset,
            label_delimiter=label_delimiter,
            pop_names=pops,
            popmax=True,
            description_text=description_text,
        )
    )

    return vcf_info_dict


def populate_info_dict(
    info_dict: Dict[str, Dict[str, str]] = INFO_DICT,
    subset_list: List[str] = SUBSETS,
    subset_pops: Dict[str, str] = POPS,
    gnomad_pops: Dict[str, str] = POPS,
    faf_pops: Dict[str, str] = FAF_POPS,
    sexes: List[str] = SEXES,
    in_silico_dict: Dict[str, Dict[str, str]] = IN_SILICO_ANNOTATIONS_INFO_DICT,
    label_delimiter: str = '_',
    bin_edges: Dict[str, str] = None,
    age_hist_data: str = None,
) -> Dict[str, Dict[str, str]]:
    """
    Call `make_info_dict` and `make_hist_dict` to populate INFO dictionary with specific sexes, population names,
    and filtering allele frequency (faf) pops.

    Used during VCF export.

    Creates:
        - INFO fields for age histograms (bin freq, n_smaller, and n_larger for heterozygous and homozygous variant carriers)
        - INFO fields for popmax AC, AN, AF, nhomalt, and popmax population
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample population, sex both for adj and raw data
        - INFO fields for filtering allele frequency (faf) annotations
        - INFO fields for variant histograms (hist_bin_freq for each histogram and hist_n_larger for DP histograms)

    :param bin_edges: Dictionary of variant annotation histograms and their associated bin edges.
    :param age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :param info_dict: INFO dict to be populated.
    :param subset_list: List of sample subsets in dataset. Default is SUBSETS.
    :param subset_pops: Dict of sample global population names to use for all subsets in `subset_list` unless the subset
        is 'gnomad', in that case `gnomad_pops` is used. Default is POPS.
    :param gnomad_pops: Dict of sample global population names for gnomAD genomes. Default is POPS.
    :param faf_pops: Dict with faf pop names (keys) and descriptions (values).  Default is FAF_POPS.
    :param sexes: gnomAD sample sexes used in VCF export. Default is SEXES.
    :param in_silico_dict: Dictionary of in silico predictor score descriptions.
    :param label_delimiter: String to use as delimiter when making group label combinations.
    :return: Updated INFO dictionary for VCF export.
    """
    vcf_info_dict = info_dict.copy()

    # Remove MISSING_INFO_FIELDS from info dict
    for field in MISSING_INFO_FIELDS:
        vcf_info_dict.pop(field, None)

    # Add allele-specific fields to info dict, including AS_VQSR_FIELDS
    vcf_info_dict.update(
        add_as_info_dict(info_dict=info_dict, as_fields=AS_FIELDS + AS_VQSR_FIELDS)
    )

    for subset in subset_list:
        if subset == 'gnomad':
            description_text = ' in gnomAD'
            pops = gnomad_pops
        else:
            description_text = '' if subset == '' else f' in {subset} subset'
            pops = subset_pops

        vcf_info_dict.update(
            populate_subset_info_dict(
                subset=subset,
                description_text=description_text,
                pops=pops,
                faf_pops=faf_pops,
                sexes=sexes,
                label_delimiter=label_delimiter,
            )
        )

    if age_hist_data:
        age_hist_data = '|'.join(str(x) for x in age_hist_data)

    vcf_info_dict.update(
        make_info_dict(
            prefix='',
            label_delimiter=label_delimiter,
            bin_edges=bin_edges,
            popmax=True,
            age_hist_data=age_hist_data,
        )
    )

    # Add in silico prediction annotations to info_dict
    vcf_info_dict.update(in_silico_dict)

    return vcf_info_dict


def make_info_expr(
    t: Union[hl.MatrixTable, hl.Table], hist_prefix: str = '',
) -> Dict[str, hl.expr.Expression]:
    """
    Make Hail expression for variant annotations to be included in VCF INFO field.

    :param t: Table/MatrixTable containing variant annotations to be reformatted for VCF export.
    :param hist_prefix: Prefix to use for histograms.
    :return: Dictionary containing Hail expressions for relevant INFO annotations.
    :rtype: Dict[str, hl.expr.Expression]
    """
    vcf_info_dict = {}
    # Add site-level annotations to vcf_info_dict
    for field in SITE_FIELDS:
        vcf_info_dict[field] = t['release_ht_info'][f'{field}']

    # Add AS annotations to info dict
    for field in AS_FIELDS:
        vcf_info_dict[field] = t['release_ht_info'][f'{field}']
    for field in VQSR_FIELDS:
        vcf_info_dict[field] = t['vqsr'][f'{field}']

    # Add region_flag and allele_info fields to info dict
    for field in ALLELE_TYPE_FIELDS:
        vcf_info_dict[field] = t['allele_info'][f'{field}']
    for field in REGION_FLAG_FIELDS:
        vcf_info_dict[field] = t['region_flag'][f'{field}']

    # Add underscore to hist_prefix if it isn't empty
    if hist_prefix != '':
        hist_prefix += '_'

    # Histograms to export are:
    # gq_hist_alt, gq_hist_all, dp_hist_alt, dp_hist_all, ab_hist_alt
    # We previously dropped:
    # _n_smaller for all hists
    # _bin_edges for all hists
    # _n_larger for all hists EXCEPT DP hists
    for hist in HISTS:
        hist_type = f'{hist_prefix}qual_hists'
        hist_dict = {
            f'{hist}_bin_freq': hl.delimit(t[hist_type][hist].bin_freq, delimiter='|'),
        }
        vcf_info_dict.update(hist_dict)

        if 'dp' in hist:
            vcf_info_dict.update({f'{hist}_n_larger': t[hist_type][hist].n_larger},)

    # Add in silico annotations to info dict
    vcf_info_dict['cadd_raw_score'] = t['cadd']['raw_score']
    vcf_info_dict['cadd_phred'] = t['cadd']['phred']

    vcf_info_dict['revel_score'] = t['revel']['revel_score']

    # In the v3.1 release in silico files this was max_ds, but changed to splice_ai_score in releases after v3.1
    vcf_info_dict['splice_ai_max_ds'] = t['splice_ai']['splice_ai_score']
    vcf_info_dict['splice_ai_consequence'] = t['splice_ai']['splice_consequence']

    vcf_info_dict['primate_ai_score'] = t['primate_ai']['primate_ai_score']

    return vcf_info_dict


def unfurl_nested_annotations(
    t: Union[hl.MatrixTable, hl.Table],
    full_release: bool = True,
    subset_release: bool = False,
    full_for_subset: bool = False,
) -> hl.expr.StructExpression:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from 
    variant annotation arrays, where the
    values of the dictionary are Hail Expressions describing how to access the 
    corresponding values.

    .. note::

       One and only one of `full_release`, `subset_release`, 
       or `full_for_subset` must be True.

    If `full_release` is True the following will be unfurled:
        - frequencies
        - popmax
        - faf
        - age histograms

    If `subset_release` is True the following will be unfurled (`hgdp_tgp_freq` 
    prefix on the freq annotation):
       - frequencies

    If `full_for_subset` is True the following will be unfurled (expects 
    'full' prefix on these annotations):
        - frequencies
        - popmax
        - faf

    :param t: Table/MatrixTable containing the nested variant annotation arrays to be 
    unfurled.
    :param full_release: Whether to unfurl gnomAD frequencies, popmax, faf, and age 
    histograms for the full release. Default is True.
    :param subset_release: Whether to unfurl frequencies for a subset. Default 
    is False.
    :param full_for_subset: Whether to unfurl full release frequencies, 
    popmax, and faf for addition to a subset release. Default is False.
    :param entries_to_remove: Optional Set of frequency entries to remove for 
    vcf_export.
    :return: StructExpression containing variant annotations and their corresponding 
    expressions and updated entries and set of frequency entries to remove
        to remove from the VCF.
    """
    expr_dict = {}

    if (full_release + subset_release + full_for_subset) != 1:
        raise ValueError(
            'One and only one of full_release, subset_release,'
            ' or full_for_subset must be set to True'
        )

    prefix = ''

    # Setting prefix with '_' as delimiter for obtaining globals
    if full_for_subset:
        prefix = 'full_'
    if subset_release:
        prefix = 'subset_'

    popmax = f'{prefix}popmax'
    faf = f'{prefix}faf'
    faf_idx = hl.eval(t.globals[f'{prefix}faf_index_dict'])

    # Unfurl freq index dict
    # Cycles through each key and index (e.g., k=adj_afr, i=31)
    logger.info('Unfurling freq data...')

    # Resetting prefix with '-' as delimiter to match values in freq_idx and faf_idx
    if full_for_subset:
        prefix = f'full-'

    if full_release or full_for_subset:
        logger.info('Adding popmax data...')
        combo_dict = {
            f'{prefix}popmax': t[popmax].pop,
            f'{prefix}AC-popmax': t[popmax].AC,
            f'{prefix}AN-popmax': t[popmax].AN,
            f'{prefix}AF-popmax': t[popmax].AF,
            f'{prefix}nhomalt-popmax': t[popmax].homozygote_count,
            f'{prefix}faf95-popmax': t[popmax].faf95,
        }
        expr_dict.update(combo_dict)

        logger.info('Unfurling faf data...')
        for (
            k,
            i,
        ) in (
            faf_idx.items()
        ):  # NOTE: faf annotations are all done on adj-only groupings
            combo_dict = {
                f'{prefix}faf95-{k}': t[faf][i].faf95,
                f'{prefix}faf99-{k}': t[faf][i].faf99,
            }
            expr_dict.update(combo_dict)

    return hl.struct(**expr_dict)


def filter_to_test(
    t: Union[hl.Table, hl.MatrixTable], num_partitions: int = 2
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter Table/MatrixTable to `num_partitions` partitions on chr20, chrX, and chrY for testing.

    :param t: Input Table/MatrixTable to filter.
    :param num_partitions: Number of partitions to grab from each chromosome.
    :return: Input Table/MatrixTable filtered to `num_partitions` on chr20, chrX, and chrY.
    """
    logger.info(
        'Filtering to %d partitions on chr20, chrX, and chrY (for tests only)...',
        num_partitions,
    )
    t_chr20 = hl.filter_intervals(t, [hl.parse_locus_interval('chr20')])
    t_chr20 = t_chr20._filter_partitions(range(num_partitions))

    t_chrx = hl.filter_intervals(t, [hl.parse_locus_interval('chrX')])
    t_chrx = t_chrx._filter_partitions(range(num_partitions))

    t_chry = hl.filter_intervals(t, [hl.parse_locus_interval('chrY')])
    t_chry = t_chry._filter_partitions(range(num_partitions))

    if isinstance(t, hl.MatrixTable):
        return (
            t_chr20.union_rows(t_chrx, t_chry)
            if isinstance(t, hl.MatrixTable)
            else t_chr20.union(t_chrx, t_chry)
        )
    else:
        t = t_chr20.union(t_chrx, t_chry)

    return t


def _prepare_vcf_ht(
    ht: hl.Table,
    is_subset: bool,
    freq_entries_to_remove: Set[str] = None,
    vcf_info_reorder: Optional[List[str]] = None,
) -> hl.Table:
    """
    Prepare the Table used for validity checks and VCF export.

    :param ht: Table containing the nested variant annotation arrays to be unfurled.
    :param is_subset: Whether this is for the release of a subset.
    :param freq_entries_to_remove: Frequency entries to remove for vcf_export.
    :param vcf_info_reorder: Optional list of INFO fields to reorder, the rest of the fields are added after this list.
    :return: Prepared HT for validity checks and VCF export
    """
    logger.info('Starting preparation of VCF HT...')
    logger.info('Adding non-PAR annotation...')

    if not is_subset:
        logger.info('Unfurling full nested frequency annotations and add to INFO field...')
        info_struct = unfurl_nested_annotations(ht)
    else:
        logger.info('Unfurling nested subset frequency annotations and add to INFO field...')
        info_struct_subset = unfurl_nested_annotations(
            ht,
            full_release=False,
            subset_release=True,
        )

        logger.info('Adding full gnomAD callset frequency annotations to INFO field...')
        info_struct_full_for_subset = unfurl_nested_annotations(
            ht,
            full_release=False,
            full_for_subset=True,
        )
        info_struct = info_struct_subset.annotate(**info_struct_full_for_subset)

    # NOTE: Merging rsid set into a semi-colon delimited string
    # dbsnp might have multiple identifiers for one variant
    # thus, rsid is a set annotation, starting with version b154 for dbsnp resource:
    # https://github.com/broadinstitute/gnomad_methods/blob/master/gnomad/resources/grch38/reference_data.py#L136
    # `export_vcf` expects this field to be a string, and vcf specs
    # say this field may be delimited by a semi-colon:
    # https://samtools.github.io/hts-specs/VCFv4.2.pdf
    # The v3.1 release chose only one of the rsids to keep so this also handles the case where rsid is a str.
    # Releases after v3.1 use the set format.
    logger.info('Reformatting rsid...')
    if isinstance(ht.rsid, hl.expr.SetExpression):
        rsid_expr = hl.str(';').join(ht.rsid)
    else:
        rsid_expr = ht.rsid

    logger.info('Constructing INFO field')
    ht = ht.annotate(
        region_flag=region_flag_expr(
            ht,
            prob_regions={
                'lcr': hl.read_table(resources.LCR_INTERVALS_HT), 
                'segdup': hl.read_table(resources.SEG_DUP_INTERVALS_HT)
            },
        ),
        release_ht_info=ht.info,
        info=info_struct,
        rsid=rsid_expr,
    )

    if freq_entries_to_remove:
        ht = ht.annotate_globals(
            freq_entries_to_remove=freq_entries_to_remove
        )
    else:
        ht = ht.annotate_globals(
            freq_entries_to_remove=hl.empty_set(hl.tstr),
        )

    # Select relevant fields for VCF export
    ht = ht.select('info', 'filters', 'rsid')

    if vcf_info_reorder:
        logger.info('Rearranging fields to desired order...')
        ht = ht.annotate(
            info=ht.info.select(*vcf_info_reorder, *ht.info.drop(*vcf_info_reorder))
        )

    return ht


def prepare_vcf_header_dict(
    t: Union[hl.Table, hl.MatrixTable],
    subset_list: List[str],
    pops: Dict[str, str],
    filtering_model_field: str = 'filtering_model',
    format_dict: Dict[str, Dict[str, str]] = FORMAT_DICT,
    inbreeding_coeff_cutoff: float = INBREEDING_COEFF_HARD_CUTOFF,
) -> Dict[str, Dict[str, str]]:
    """
    Prepare VCF header dictionary.

    :param t: Input MatrixTable/Table
    :param age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :param subset_list: List of sample subsets in dataset.
    :param pops: List of sample global population names for gnomAD genomes.
    :param filtering_model_field: String indicating the filtering model global annotation.
    :param format_dict: Dictionary describing MatrixTable entries. Used in header for VCF export.
    :param inbreeding_coeff_cutoff: InbreedingCoeff hard filter used for variants.
    :return: Prepared VCF header dictionary.
    """
    logger.info('Making FILTER dict for VCF...')
    filter_dict = make_vcf_filter_dict(
        hl.eval(t[filtering_model_field].snv_cutoff.min_score),
        hl.eval(t[filtering_model_field].indel_cutoff.min_score),
        inbreeding_cutoff=inbreeding_coeff_cutoff,
        variant_qc_filter='AS_VQSR',
    )

    logger.info('Making INFO dict for VCF...')
    vcf_info_dict = populate_info_dict(
        subset_list=subset_list,
        subset_pops=pops,
    )

    # Adjust keys to remove adj tags before exporting to VCF
    # VCF 4.3 specs do not allow hyphens in info fields
    new_vcf_info_dict = {
        i.replace('_adj', '').replace('-', '_'): j for i, j in vcf_info_dict.items()
    }

    header_dict = {
        'info': new_vcf_info_dict,
        'filter': filter_dict,
        'format': format_dict,
    }

    return header_dict


def _cleanup_ht_for_vcf_export(ht: hl.Table) -> hl.Table:
    """
    Clean up the Table returned by `prepare_vcf_ht` so it is ready to export to VCF with `hl.export_vcf`.

    Specifically:
    - Adjust types that are incompatible with VCFs using `adjust_vcf_incompatible_types`.

    :param ht: Table returned by `prepare_vcf_ht`.
    :return: Table ready for export to VCF and a list of fixed row annotations needed for the VCF header check.
    """
    logger.info('Adjusting VCF incompatible types...')

    # The Table is already split so there are no annotations that need to be pipe delimited
    ht = adjust_vcf_incompatible_types(ht, pipe_delimited_annotations=[])

    return ht


def _build_parameter_dict(
    ht: hl.Table,
    is_public_subset: bool = False,
) -> Dict[str, Union[bool, str, List, Dict, Set, hl.Table, None]]:
    """
    Build a dictionary of parameters to export.

    Parameters differ from subset releases (e.g., HGDP + TGP) vs full release.

    :param is_public_subset: Build the parameter list specific to the HGDP + TGP 
    subset release.
    :param is_test: Uses a checkpoint path for the prepared VCF Table that adds the 
    string 'test' to the checkpoint path.
    :return: Dictionary containing parameters needed to make the release VCF.
    """
    if is_public_subset:
        parameter_dict = {
            'pops': HGDP_TGP_POPS,
            'subsets': ['', 'gnomad'],
            'drop_hists': None,
            'include_age_hists': False,
            'sample_sum_sets_and_pops': {'gnomad': POPS, '': HGDP_TGP_POPS},
            'freq_entries_to_remove': set(),
            'age_hist_data': None,
            'filtering_model_field': 'variant_filtering_model',
        }

    else:
        parameter_dict = {
            'pops': POPS,
            'subsets': SUBSET_LIST_FOR_VCF,
            'drop_hists': ['age_hist_het_bin_edges', 'age_hist_hom_bin_edges'],
            'include_age_hists': True,
            'sample_sum_sets_and_pops': {'hgdp': HGDP_POPS, 'tgp': TGP_POPS},
        }
        # Downsampling and subset entries to remove from VCF's freq export
        # Note: Need to extract the non-standard downsamplings from the freq_meta struct to the FREQ_ENTRIES_TO_REMOVE
        freq_entries_to_remove = {
            str(x['downsampling'])
            for x in hl.eval(ht.freq_meta)
            if 'downsampling' in x
        }
        freq_entries_to_remove.update(set(COHORTS_WITH_POP_STORED_AS_SUBPOP))
        parameter_dict['freq_entries_to_remove'] = freq_entries_to_remove
        parameter_dict['filtering_model_field'] = 'filtering_model'

    return parameter_dict


if __name__ == '__main__':
    main()  # pylint: disable=E1120