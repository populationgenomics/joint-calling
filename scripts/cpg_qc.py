import os
import time
from typing import List, Optional
import click
import logging
import hail as hl
from gnomad.utils.file_utils import file_exists

from qc import resources, sample_qc as qc


logger = logging.getLogger("cpg_qc")


DEFAULT_REF = 'GRCh38'
TIMESTAMP = time.strftime('%Y%m%d-%H%M')


@click.group()
@click.option('--mt', 'mt_path', required=True,
              help='path to the matrix table .mt directory (can be on gs://)')
@click.option('--bucket', 'work_bucket', required=True,
              help='path to folder for intermediate output (can be on gs://)')
@click.option('--log-dir', 'log_dirpath', required=True,
              help='local directory store Hail logs (must be local)')
@click.option('--overwrite', 'overwrite', is_flag=True)
@click.pass_context
def cli(ctx, mt_path: str, work_bucket: str, log_dirpath: str, overwrite: bool):
    ctx.ensure_object(dict)
    ctx.obj['mt_path'] = mt_path
    ctx.obj['work_bucket'] = work_bucket
    ctx.obj['log_dirpath'] = log_dirpath
    ctx.obj['overwrite'] = overwrite


@cli.command()
@click.pass_context
@click.option('--sample-map', 'sample_map_fpath',
help=\
'sample_map: path to the sample map (must be filesystem local).'
'The sample map should be tab separated with two columns.'
'The first column is the sample ID, and the second column'
'is the gVCF path. WARNING: the sample names in the gVCFs will be'
'overwritten')
def combine_gvcfs(ctx, sample_map_fpath: str):
    mt_path = ctx.obj['mt_path']
    work_bucket = ctx.obj['work_bucket']
    log_dirpath = ctx.obj['log_dirpath']

    if file_exists(mt_path):
        print(f'{mt_path} exists, reusing')
        return
    hl.init(default_reference=DEFAULT_REF,
            log=os.path.join(log_dirpath, f'gvcfs-combine-{TIMESTAMP}.log'))

    with open(sample_map_fpath) as f:
        names_and_paths = [tuple(l.strip().split('\t')) for l in f]

    hl.experimental.run_combiner(
        [p for n, p in names_and_paths],
        out_file=mt_path,
        reference_genome=DEFAULT_REF,
        use_genome_default_intervals=True,
        tmp_path=os.path.join(work_bucket, 'tmp'),
        key_by_locus_and_alleles=True,
        overwrite=True
    )


@cli.command()
@click.pass_context
def sample_qc(ctx):
    from hail.experimental.vcf_combiner.sparse_mt_utils import lgt_to_gt

    mt_path = ctx.obj['mt_path']
    work_bucket = ctx.obj['work_bucket']
    log_dirpath = ctx.obj['log_dirpath']
    overwrite = ctx.obj['overwrite']
    resources.BUCKET = work_bucket

    mt_out_path = os.path.basename(mt_path).replace('.mt', '').replace('/', '') + '.qc.mt'

    if not overwrite and file_exists(mt_out_path):
        print(f'{mt_out_path} exists, reusing')
        return

    hl.init(default_reference=DEFAULT_REF,
            log=os.path.join(log_dirpath, f'sample-qc-{TIMESTAMP}.log'))

    mt = hl.read_matrix_table(mt_path)

    ht = qc.compute_sample_qc(mt, os.path.join(work_bucket, 'checkpoints', 'qc'))
    ht.write(resources.get_sample_qc_url(work_bucket).path, overwrite=True)

    qc.compute_sex(mt).write(resources.sex, overwrite=True)


if __name__ == '__main__':
    cli()


