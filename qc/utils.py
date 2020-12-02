import logging
import os
import subprocess
# from google.cloud import secretmanager
# from variant_qc import _version
# from variant_qc.utils import run
from gnomad.resources import VersionedTableResource, TableResource
from google.cloud import storage


def file_exists(path: str):
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        path = path.replace('gs://', '').split('/', maxsplit=1)[1]
        gs = storage.Client()
        return gs.get_bucket(bucket).get_blob(path)
    return os.path.exists(path)


def gs_download_file(path):
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        path = path.replace('gs://', '').split('/', maxsplit=1)[1]
        gs = storage.Client()
        blob = gs.get_bucket(bucket).get_blob(path)
        if blob:
            return blob.download_as_string()
    else:
        if os.path.isfile(path):
            with open(path) as f:
                return f.read()
    return None


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
