#!/usr/bin/env python

"""
Generates annotations for random_forest.py
"""

import logging
from os.path import join, basename
from typing import Dict, Optional
import subprocess
import click
import hail as hl

from gnomad.sample_qc.relatedness import generate_trio_stats_expr
from gnomad.utils.annotations import (
    add_variant_type,
    annotate_adj,
)
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.sparse_mt import (
    split_info_annotation,
    split_lowqual_annotation,
)
from gnomad.utils.vep import vep_or_lookup_vep

from joint_calling.utils import file_exists
from joint_calling import utils, _version
import joint_calling

logger = logging.getLogger('qc-annotations')
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--out-allele-data-ht',
    'out_allele_data_ht_path',
    required=True,
)
@click.option(
    '--out-qc-ac-ht',
    'out_qc_ac_ht_path',
    required=True,
)
@click.option(
    '--out-fam-stats-ht',
    'out_fam_stats_ht_path',
    required=True,
)
@click.option(
    '--out-vep-ht',
    'out_vep_ht_path',
)
@click.option(
    '--mt',
    'mt_path',
    required=True,
    callback=utils.get_validation_callback(ext='mt', must_exist=True),
    help='path to the input MatrixTable',
)
@click.option(
    '--hard-filtered-samples-ht',
    'hard_filtered_samples_ht_path',
    required=True,
    help='Path to a table with only samples that passed filters '
    '(it\'s generated by sample QC)',
)
@click.option('--fam-file', 'trios_fam_ped_file', help='PED file for the samples')
@click.option(
    '--meta-ht',
    'meta_ht_path',
    required=True,
    help='',
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
    '--vep-version',
    'vep_version',
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,too-many-statements,missing-function-docstring
    out_allele_data_ht_path: str,
    out_qc_ac_ht_path: str,
    out_fam_stats_ht_path: str,
    out_vep_ht_path: Optional[str],
    mt_path: str,
    hard_filtered_samples_ht_path: str,
    trios_fam_ped_file: Optional[str],
    meta_ht_path: str,
    work_bucket: str,
    local_tmp_dir: str,
    overwrite: bool,
    vep_version: Optional[str],
):
    utils.init_hail('qc_annotations', local_tmp_dir)

    all_samples_mt = utils.get_mt(mt_path)
    hard_filtered_mt = utils.get_mt(
        mt_path,
        hard_filtered_samples_to_remove_ht=hl.read_table(hard_filtered_samples_ht_path),
        meta_ht=hl.read_table(meta_ht_path),
        add_meta=True,
    )

    generate_allele_data(
        ht=hard_filtered_mt.rows(),
        out_ht_path=out_allele_data_ht_path,
        overwrite=overwrite,
    )

    qc_ac_ht = generate_ac(
        mt=hard_filtered_mt,
        out_ht_path=out_qc_ac_ht_path,
        overwrite=overwrite,
    )

    if not trios_fam_ped_file:
        trios_fam_ped_file = _make_fam_file(
            sex_ht=hl.read_table(meta_ht_path),
            work_bucket=work_bucket,
        )

    fam_stats_ht = generate_fam_stats(
        hard_filtered_mt,
        out_fam_stats_ht_path,
        overwrite=overwrite,
        trios_fam_ped_file=trios_fam_ped_file,
    )

    if fam_stats_ht:
        export_transmitted_singletons_vcf(
            fam_stats_ht=fam_stats_ht,
            qc_ac_ht=qc_ac_ht,
            work_bucket=work_bucket,
            overwrite=overwrite,
        )

    if vep_version or out_vep_ht_path:
        if not out_vep_ht_path:
            logger.critical('--out-vep-ht must be specified along with --vep-version')
        run_vep(
            out_ht_path=str(out_vep_ht_path),
            mt=all_samples_mt,
            vep_version=vep_version,
            work_bucket=work_bucket,
            overwrite=overwrite,
        )


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


def generate_allele_data(
    ht: hl.Table,
    out_ht_path: str,
    overwrite: bool,
) -> hl.Table:
    """
    Returns bi-allelic sites HT with the following annotations:
     - allele_data (nonsplit_alleles, has_star, variant_type, and n_alt_alleles)
    :param Table ht: full unsplit HT
    :param Table out_ht_path: path to writ the resulting table
    :param overwrite: overwrite checkpoints if they exist
    :return: Table with allele data annotations
    """
    logger.info('Generate allele data')
    if not overwrite and file_exists(out_ht_path):
        return hl.read_table(out_ht_path)

    ht = ht.select()
    allele_data = hl.struct(
        nonsplit_alleles=ht.alleles, has_star=hl.any(lambda a: a == '*', ht.alleles)
    )
    ht = ht.annotate(allele_data=allele_data.annotate(**add_variant_type(ht.alleles)))

    ht = hl.split_multi_hts(ht)
    ht = ht.filter(hl.len(ht.alleles) > 1)
    allele_type = (
        hl.case()
        .when(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv')
        .when(hl.is_insertion(ht.alleles[0], ht.alleles[1]), 'ins')
        .when(hl.is_deletion(ht.alleles[0], ht.alleles[1]), 'del')
        .default('complex')
    )
    ht = ht.annotate(
        allele_data=ht.allele_data.annotate(
            allele_type=allele_type, was_mixed=ht.allele_data.variant_type == 'mixed'
        )
    )
    ht.write(out_ht_path, overwrite=True)
    return ht


def generate_ac(mt: hl.MatrixTable, out_ht_path: str, overwrite: bool) -> hl.Table:
    """
    Creates Table containing allele counts per variant.
    :param Table mt: input hard-filtered, meta-annotated matrix table
    :param Table out_ht_path: path to writ the resulting table
    :param overwrite: overwrite checkpoints if they exist
    :return: table containing allele counts annotations:
        - `ac_qc_samples_raw`: Allele count of high quality samples
        - `ac_qc_samples_unrelated_raw`: Allele count of high quality
           unrelated samples
        - `ac_release_samples_raw`: Allele count of release samples
        - `ac_qc_samples_adj`: Allele count of high quality samples after
           adj filtering
        - `ac_qc_samples_unrelated_adj`: Allele count of high quality
           unrelated samples after adj filtering
        - `ac_release_samples_adj`: Allele count of release samples after
           adj filtering
    """
    logger.info('Generate AC per variant')
    if not overwrite and file_exists(out_ht_path):
        logger.info(f'Reusing {out_ht_path}')
        return hl.read_table(out_ht_path)

    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
    mt = mt.filter_cols(mt.meta.high_quality)
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    mt = annotate_adj(mt)
    mt = mt.annotate_rows(
        ac_qc_samples_raw=hl.agg.sum(mt.GT.n_alt_alleles()),
        ac_qc_samples_unrelated_raw=hl.agg.filter(
            ~mt.meta.related,
            hl.agg.sum(mt.GT.n_alt_alleles()),
        ),
        ac_release_samples_raw=hl.agg.filter(
            mt.meta.release, hl.agg.sum(mt.GT.n_alt_alleles())
        ),
        ac_qc_samples_adj=hl.agg.filter(mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
        ac_qc_samples_unrelated_adj=hl.agg.filter(
            ~mt.meta.related & mt.adj,
            hl.agg.sum(mt.GT.n_alt_alleles()),
        ),
        ac_release_samples_adj=hl.agg.filter(
            mt.meta.release & mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())
        ),
    )
    ht = mt.rows()
    ht = ht.repartition(10000, shuffle=False)
    ht.write(out_ht_path, overwrite=True)
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


def generate_fam_stats(
    mt: hl.MatrixTable,
    out_fam_stats_ht_path: str,
    overwrite: bool,
    trios_fam_ped_file: str,
) -> hl.Table:
    """
    Calculate transmission and de novo mutation statistics using trios in the dataset.
    :param mt: input MatrixTable
    :param out_fam_stats_ht_path: path to write the resulting Table to
    :param work_bucket: bucket to write intermediate and output files to
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
        f'Generating family stats using {mt.count_cols()} samples from {len(ped.trios)} trios.'
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


def export_transmitted_singletons_vcf(
    fam_stats_ht: hl.Table,
    qc_ac_ht: hl.Table,
    work_bucket: str,
    overwrite: bool = False,
) -> Dict[str, str]:
    """
    Exports the transmitted singleton Table to a VCF.
    :return: None
    """
    output_vcf_paths = {
        conf: join(work_bucket, f'transmitted-singletons-{conf}.vcf.bgz')
        for conf in ['adj', 'raw']
    }
    if not overwrite and all(file_exists(path) for path in output_vcf_paths.values()):
        return output_vcf_paths

    for transmission_confidence in ['raw', 'adj']:
        ts_ht = qc_ac_ht.filter(
            (
                fam_stats_ht[qc_ac_ht.key][f'n_transmitted_{transmission_confidence}']
                == 1
            )
            & (qc_ac_ht.ac_qc_samples_raw == 2)
        )
        ts_ht = ts_ht.annotate(s=hl.null(hl.tstr))
        ts_mt = ts_ht.to_matrix_table_row_major(columns=['s'], entry_field_name='s')
        ts_mt = ts_mt.filter_cols(False)
        hl.export_vcf(
            ts_mt,
            output_vcf_paths[transmission_confidence],
            tabix=True,
        )
    return output_vcf_paths


def run_vep(
    out_ht_path: str,
    mt: hl.MatrixTable,
    vep_version: Optional[str],
    work_bucket: str,
    overwrite: bool = False,
) -> hl.Table:
    """
    Returns a table with a VEP annotation for each variant in the raw MatrixTable.
    :param mt: keyed by locus and allele, soft-filtered
    :param vep_version:
    :return: VEPed Table
    """
    if not overwrite and file_exists(out_ht_path):
        return hl.read_table(out_ht_path)
    vep_config_local = join(joint_calling.get_package_path(), 'vep-config.json')
    vep_config_gs = join(work_bucket, basename(vep_config_local))
    subprocess.run(
        f'gsutil cp {vep_config_local} {vep_config_gs}', check=False, shell=True
    )
    ht = mt.rows()
    ht = ht.filter(hl.len(ht.alleles) > 1)
    ht = hl.split_multi_hts(ht)
    # TODO: rerun VEP instead of reading gnomAD vep reference file
    ht = vep_or_lookup_vep(ht, vep_version=vep_version, vep_config_path=vep_config_gs)
    ht = ht.annotate_globals(version=f'v{vep_version}')
    ht.write(out_ht_path, overwrite=True)
    return ht


if __name__ == '__main__':
    main()  # pylint: disable=E1120
