import csv
import os
import time
from typing import List, Optional
import click
import logging
import pandas as pd
import hail as hl
from gnomad.utils.file_utils import file_exists

from qc import resources, sample_qc as qc
from qc.utils import safe_mkdir

logger = logging.getLogger("cpg_qc")

DEFAULT_REF = 'GRCh38'
TIMESTAMP = time.strftime('%Y%m%d-%H%M')


def init_hail(name, local_tmp_dir):
    hl_log = os.path.join(safe_mkdir(os.path.join(local_tmp_dir, 'log')), f'{name}-{TIMESTAMP}.log')
    hl.init(default_reference=DEFAULT_REF, log=hl_log)


@click.group()
@click.option('--mt', 'mt_path', required=True,
              help='path to the matrix table .mt directory (can be on gs://)')
@click.option('--bucket', 'work_bucket', required=True,
              help='path to folder for intermediate output (can be on gs://)')
@click.option('--local-tmp-dir', 'local_tmp_dir', required=True,
              help='local directory to store temporary files and Hail logs (must be local)')
@click.option('--reuse', 'reuse', is_flag=True)
@click.pass_context
def cli(ctx: click.core.Context,
        mt_path: str,
        work_bucket: str,
        local_tmp_dir: str,
        reuse: bool):
    ctx.ensure_object(dict)
    ctx.obj['mt_path'] = mt_path
    ctx.obj['work_bucket'] = work_bucket
    ctx.obj['local_tmp_dir'] = local_tmp_dir
    ctx.obj['reuse'] = reuse
    resources.BUCKET = work_bucket


@cli.command()
@click.pass_context
@click.option('--sample-map', 'sample_map_csv',
help=\
'sample_map: path to the sample map (must be filesystem local).'
'The sample map should be comma-separated with the following columns:\n'
'sample,gvcfs,contamination,alignment_summary_metrics,duplicate_metrics,insert_size_metrics,wgs_metrics\n'
'The first column is the sample ID.')
def combine_gvcfs(ctx: click.core.Context, sample_map_csv: str):
    mt_path = ctx.obj['mt_path']
    work_bucket = ctx.obj['work_bucket']
    local_tmp_dir = ctx.obj['local_tmp_dir']
    reuse = ctx.obj['reuse']

    if reuse and file_exists(mt_path):
        logger.info(f'Reusing existing output: {mt_path}')
        return

    try:
        init_hail('combine_gvcfs', local_tmp_dir)
    except:
        pass

    sample_df = pd.read_csv(sample_map_csv, sep=',')

    hl.experimental.run_combiner(
        sample_df['gvcf'],
        out_file=mt_path,
        reference_genome=DEFAULT_REF,
        use_genome_default_intervals=True,
        tmp_path=os.path.join(work_bucket, 'tmp'),
        key_by_locus_and_alleles=True,
        overwrite=True
    )


@cli.command()
@click.pass_context
@click.option('--out-ht', 'out_ht_path', required=True)
def sample_qc(ctx: click.core.Context, out_ht_path: str):
    mt_path = ctx.obj['mt_path']
    work_bucket = ctx.obj['work_bucket']
    local_tmp_dir = ctx.obj['local_tmp_dir']
    reuse = ctx.obj['reuse']

    if reuse and file_exists(out_ht_path):
        logger.info(f'Reusing existing output: {out_ht_path}')
        return

    init_hail('sample_qc', local_tmp_dir)

    mt = hl.read_matrix_table(mt_path)
    basename = os.path.splitext(out_ht_path)[0]
    ht = qc.compute_sample_qc(mt, basename)
    ht.write(out_ht_path, overwrite=True)


@cli.command()
@click.pass_context
@click.option('--out-mt', 'out_mt_path', required=True)
def compute_qc_mt(ctx: click.core.Context, out_mt_path: str):
    mt_path = ctx.obj['mt_path']
    work_bucket = ctx.obj['work_bucket']
    local_tmp_dir = ctx.obj['local_tmp_dir']
    reuse = ctx.obj['reuse']
    if reuse and file_exists(out_mt_path):
        logger.info(f'Reusing existing output: {out_mt_path}')
        return

    init_hail('compute_qc_mt', local_tmp_dir)
    mt = hl.read_matrix_table(mt_path)

    mt = qc.compute_qc_mt(mt)
    mt.write(out_mt_path, overwrite=True)


@cli.command()
@click.pass_context
@click.option('--out-ht', 'out_ht_path', required=True)
def impute_sex(ctx: click.core.Context, out_ht_path: str):
    mt_path = ctx.obj['mt_path']
    work_bucket = ctx.obj['work_bucket']
    local_tmp_dir = ctx.obj['local_tmp_dir']
    reuse = ctx.obj['reuse']

    if reuse and file_exists(out_ht_path):
        logger.info(f'Reusing existing output: {out_ht_path}')
        return

    init_hail('impute_sex', local_tmp_dir)

    mt = hl.read_matrix_table(mt_path)

    from gnomad.resources.grch38 import telomeres_and_centromeres
    ht = qc.annotate_sex(
        mt,
        excluded_intervals=telomeres_and_centromeres.ht(),
        aaf_threshold=0.001,
        f_stat_cutoff=0.5,
        gt_expr="LGT",
    )
    ht.write(out_ht_path, overwrite=True)


@cli.command()
@click.pass_context
@click.option('--out-ht', 'out_ht_path', required=True)
@click.option('--sex-ht', 'sex_ht_path', required=True)
@click.option('--biallelic-qc-ht', 'biallelic_qc_ht_path', required=True)
@click.option('--min-cov', 'min_cov', default=18)
@click.option('--sample-map', 'sample_map_csv',
help=\
'sample_map: path to the sample map (must be filesystem local).'
'The sample map should be comma-separated with the following columns:\n'
'sample,gvcfs,contamination,alignment_summary_metrics,duplicate_metrics,insert_size_metrics,wgs_metrics\n'
'The first column is the sample ID.')
def compute_hard_filters(
        ctx: click.core.Context,
        out_ht_path: str,
        sex_ht_path: str,
        biallelic_qc_ht_path: str,
        min_cov: int,
        sample_map_csv: str
    ):
    mt_path = ctx.obj['mt_path']
    work_bucket = ctx.obj['work_bucket']
    local_tmp_dir = ctx.obj['local_tmp_dir']
    reuse = ctx.obj['reuse']

    if reuse and file_exists(out_ht_path):
        logger.info(f'Reusing existing output: {out_ht_path}')
        return

    init_hail('compute_hard_filters', local_tmp_dir)
    mt = hl.read_matrix_table(mt_path)

    sample_df = pd.read_csv(sample_map_csv, sep=',')
    metrics_ht = qc.parse_metrics(sample_df, local_tmp_dir)

    qc.compute_hard_filters(
        mt,
        cov_threshold=min_cov,
        biallelic_qc_ht=hl.read_table(biallelic_qc_ht_path),
        sex_ht=hl.read_table(sex_ht_path),
        metrics_ht=metrics_ht
    ).write(out_ht_path, overwrite=True)


@cli.command()
@click.pass_context
@click.option('--out-ht', 'out_ht_path', required=True)
def run_pc_relate(
        ctx: click.core.Context,
        out_ht_path: str,
    ):
    mt_path = ctx.obj['mt_path']
    work_bucket = ctx.obj['work_bucket']
    local_tmp_dir = ctx.obj['local_tmp_dir']
    reuse = ctx.obj['reuse']
    if reuse and file_exists(out_ht_path):
        logger.info(f'Reusing existing output: {out_ht_path}')
        return

    init_hail('run_pc_relate', local_tmp_dir)
    qc_mt = hl.read_matrix_table(mt_path)
    eig, scores, _ = hl.hwe_normalized_pca(qc_mt.LGT, k=10, compute_loadings=False)
    # scores = scores.checkpoint(v3_pc_relate_pca_scores.path, overwrite=reuse, _read_if_exists=not reuse)
    relatedness_ht = hl.pc_relate(
        qc_mt.LGT,
        min_individual_maf=0.01,
        scores_expr=scores[qc_mt.col_key].scores,
        block_size=4096,
        min_kinship=0.05,
        statistics='all')
    relatedness_ht.write(out_ht_path, overwrite=True)


if __name__ == '__main__':
    cli()


