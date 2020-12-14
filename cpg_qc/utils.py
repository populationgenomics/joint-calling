import logging
import os
import subprocess
import sys
import time
from os.path import isdir, isfile, exists
from google.cloud import storage
import hashlib


def file_exists(path: str):
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        path = path.replace('gs://', '').split('/', maxsplit=1)[1]
        if path.endswith('.mt') or path.endswith('.mt/'):
            path = os.path.join(path, '_SUCCESS')
        gs = storage.Client()
        return gs.get_bucket(bucket).get_blob(path)
    return os.path.exists(path)


def gs_cache_file(fpath: str, local_tmp_dir: str = None) -> str:
    """
    :param fpath: local or gs:// path. If the latter, the file will be downloaded
                  and cached if local_tmp_dir is provided, the local path will be
                  returned
    :param local_tmp_dir: a directory to cache files downloaded from cloud
    :return: file path
    """
    if fpath.startswith('gs://'):
        fname = os.path.basename(fpath) + '_' + hashlib.md5(fpath.encode()).hexdigest()[:6]
        local_fpath = os.path.join(local_tmp_dir, fname)
        if not exists(local_fpath):
            bucket = fpath.replace('gs://', '').split('/')[0]
            path = fpath.replace('gs://', '').split('/', maxsplit=1)[1]
            gs = storage.Client()
            blob = gs.get_bucket(bucket).get_blob(path)
            if blob:
                blob.download_to_filename(local_fpath)
    else:
        local_fpath = fpath
    return local_fpath


def safe_mkdir(dirpath, descriptive_name=''):
    """ Multiprocessing-safely and recursively creates a directory
    """
    if not dirpath:
        sys.stderr.write(f'Path is empty: {descriptive_name if descriptive_name else ""}\n')

    if isdir(dirpath):
        return dirpath

    if isfile(dirpath):
        sys.stderr.write(descriptive_name + ' ' + dirpath + ' is a file.\m')

    num_tries = 0
    max_tries = 10

    while not exists(dirpath):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dirpath)
        except OSError as e:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dirpath


def run(cmd, silent=False):
    """Run the provided command, logging details and checking for errors.
    """
    if not silent:
        logging.warning(' '.join(str(x) for x in cmd) if not isinstance(cmd, str) else cmd)
    subprocess.check_call(cmd, shell=True, executable=find_bash())


def find_bash():
    for test_bash in ["/bin/bash", "/usr/bin/bash", "/usr/local/bin/bash"]:
        if test_bash and os.path.exists(test_bash):
            return test_bash
    raise IOError("Could not find bash in any standard location. Needed for unix pipes")



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
