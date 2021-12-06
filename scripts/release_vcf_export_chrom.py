import logging
import pickle

import click
import hail as hl

from gnomad.resources.grch38.gnomad import (
    COHORTS_WITH_POP_STORED_AS_SUBPOP,
    HGDP_POPS,
    TGP_POPS,
    TGP_POP_NAMES,
    POPS,
    SUBSETS,
)
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.utils.filtering import remove_fields_from_constant
from gnomad.utils.vcf import (
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    build_vcf_export_reference,
    FAF_POPS,
    FORMAT_DICT,
    rekey_new_reference,
    REGION_FLAG_FIELDS,
    RF_FIELDS,
    SITE_FIELDS,
)
from joint_calling import utils, _version
from joint_calling.utils import get_validation_callback


logger = logging.getLogger(__file__)
logger.setLevel('INFO')


# Add new site fields
NEW_SITE_FIELDS = [
    'monoallelic',
    'transmitted_singleton',
]
SITE_FIELDS.extend(NEW_SITE_FIELDS)

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

# VCF INFO field reordering
VCF_INFO_REORDER = ['AC-adj', 'AN-adj', 'AF-adj', 'popmax', 'faf95-popmax']
HGDP_TGP_VCF_INFO_REORDER = [
    'AC-adj',
    'AN-adj',
    'AF-adj',
    'AC-raw',
    'AN-raw',
    'AF-raw',
    'gnomad-AC-adj',
    'gnomad-AN-adj',
    'gnomad-AF-adj',
    'gnomad-popmax',
    'gnomad-faf95-popmax',
    'gnomad-AC-raw',
    'gnomad-AN-raw',
    'gnomad-AF-raw',
]


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--ht',
    'ht_path',
    required=True,
    callback=get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--vcf-header-txt',
    'vcf_header_txt_path',
    required=True,
    callback=get_validation_callback(ext='txt', must_exist=True),
)
@click.option(
    '--out-vcf',
    'out_vcf_path',
    required=True,
    callback=get_validation_callback(ext='vcf.bgz'),
)
@click.option(
    '--name',
    'name',
    required=True,
    help='name of the dataset',
)
@click.option(
    '--chromosome',
    'chromosome',
    help='write a VCF for a chromosome',
)
@click.option(
    '--local-tmp-dir',
    'local_tmp_dir',
    help='local directory for temporary files and Hail logs (must be local).',
)
def main(
    ht_path: str,
    vcf_header_txt_path: str,
    out_vcf_path: str,
    name: str,
    chromosome: str,
    local_tmp_dir: str,
):  # pylint: disable=missing-function-docstring
    utils.init_hail(__file__, local_tmp_dir)

    ht = hl.read_table(ht_path)

    logger.info('Loading VCF header dict...')
    with hl.hadoop_open(vcf_header_txt_path, 'rb') as f:
        header_dict = pickle.load(f)

    if chromosome:
        logger.info(f'Exporting chromosome {chromosome}....')
        ht = hl.filter_intervals(ht, [hl.parse_locus_interval(chromosome)])

    export_reference = build_vcf_export_reference(name)

    hl.export_vcf(
        rekey_new_reference(ht, export_reference),
        out_vcf_path,
        metadata=header_dict,
        tabix=True,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
