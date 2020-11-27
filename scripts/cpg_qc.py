import os
import time
from typing import List
import click
# from google.cloud import secretmanager
# from variant_qc import _version
# from variant_qc.utils import run
from google.cloud import storage
import hail as hl


DEFAULT_REF = 'GRCh38'
TIMESTAMP = time.strftime('%Y%m%d-%H%M')


@click.group()
@click.option('--mt-file', 'mt_fpath', required=True,
              help='path to matrix table (can be on gs://)')
@click.option('--bucket', 'work_bucket', required=True,
              help='path to folder for intermediate output (can be on gs://)')
@click.option('--log-dir', 'log_dirpath', required=True,
              help='local directory store Hail logs (must be local)')
@click.pass_context
def cli(ctx, mt_fpath: str, work_bucket: str, log_dirpath: str):
    ctx.ensure_object(dict)
    ctx.obj['mt_fpath'] = mt_fpath
    ctx.obj['work_bucket'] = work_bucket
    ctx.obj['log_dirpath'] = log_dirpath


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
    mt_fpath = ctx.obj['mt_fpath']
    work_bucket = ctx.obj['work_bucket']
    log_dirpath = ctx.obj['log_dirpath']

    if file_exists(mt_fpath):
        print(f'{mt_fpath} exists, reusing')
        return
    hl.init(default_reference=DEFAULT_REF,
            log=os.path.join(log_dirpath, f'gvcfs-combine-{TIMESTAMP}.log'))

    with open(sample_map_fpath) as f:
        names_and_paths = [tuple(l.strip().split('\t')) for l in f]

    hl.experimental.run_combiner(
        [p for n, p in names_and_paths],
        out_file=mt_fpath,
        reference_genome=DEFAULT_REF,
        use_genome_default_intervals=True,
        tmp_path=os.path.join(work_bucket, 'tmp'),
        overwrite=True
    )


@cli.command()
@click.pass_context
def sample_qc(ctx):
    mt_fpath = ctx.obj['mt_fpath']
    work_bucket = ctx.obj['work_bucket']
    log_dirpath = ctx.obj['log_dirpath']

    mt_out_path = os.path.basename(mt_fpath).replace('.mt', '').replace('/', '') + '.qc.mt'

    if file_exists(mt_out_path):
        print(f'{mt_out_path} exists, reusing')
        return
    hl.init(default_reference=DEFAULT_REF,
            log=os.path.join(log_dirpath, f'sample-qc-{TIMESTAMP}.log'))

    mt = hl.read_matrix_table(mt_fpath)
    mt = hl.sample_qc(mt)
    mt.write(os.path.join(work_bucket, mt_out_path))


#     from google.cloud import storage
#     # If you don't specify credentials when constructing the client, the
#     # client library will look for credentials in the environment.
#     storage_client = storage.Client()
#     # Make an authenticated API request
#     # buckets = list(storage_client.list_buckets())
#     # print(buckets)
#     f = storage_client.get_bucket('warp-playground').get_blob('NA12878_24RG_small.hg38.bam')
#     print(f)
#
#     sys.exit(1)
#
#
#     slack_channel = _init_notifications()
#
#     # run(f'import_vcf.py '
#     #     f'--vcf {vcf} '
#     #     f'--results_bucket {results_bucket} '
#     #     f'--slack_channel {slack_channel} '
#     #     f'--overwrite '
#     #     )
#     run(f'variantqc.py '
#         f'--results_bucket {results_bucket} '
#         f'--debug '
#         f'--slack_channel {slack_channel} '
#         f'--annotate_for_rf '
#         f'--train_rf '
#         f'--apply_rf '
#         f'--finalize '
#         )
#
#
# def _init_notifications():
#     gcp_project_number = os.getenv('GCP_PROJECT_NUMBER')
#     assert gcp_project_number
#     slack_channel = os.getenv('SLACK_CHANNEL')
#     slack_token_secret_name = (
#         f'projects/{gcp_project_number}/secrets/slack_token/versions/2')
#
#     # Cache the Slack client.
#     secret_manager = secretmanager.SecretManagerServiceClient()
#     slack_token_response = secret_manager.access_secret_version(
#         request={"name": slack_token_secret_name})
#     slack_token = slack_token_response.payload.data.decode('UTF-8')
#     os.environ['SLACK_API_TOKEN'] = slack_token
#     return slack_channel


def file_exists(path: str):
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        path = path.replace('gs://', '').split('/', maxsplit=1)[1]
        gs = storage.Client()
        return gs.get_bucket(bucket).get_blob(path)
    return os.path.exists(path)


if __name__ == '__main__':
    cli()
