import os
import logging
import click
from google.cloud import secretmanager
from variant_qc import _version
from variant_qc.utils import run


@click.command()
@click.version_option(_version.__version__)
@click.argument('vcf')
@click.option('-b', '--results-bucket', 'results_bucket')
def main(vcf, results_bucket):
    slack_channel = _init_notifications()

    run(f'import_vcf.py '
        f'--vcf {vcf} '
        f'--results_bucket {results_bucket} '
        f'--slack_channel {slack_channel} '
        f'--overwrite '
    )


def _init_notifications():
    gcp_project_number = os.getenv('GCP_PROJECT_NUMBER')
    assert gcp_project_number
    slack_channel = os.getenv('SLACK_CHANNEL')
    slack_token_secret_name = (
        f'projects/{gcp_project_number}/secrets/slack_token/versions/2')

    # Cache the Slack client.
    secret_manager = secretmanager.SecretManagerServiceClient()
    slack_token_response = secret_manager.access_secret_version(
        request={"name": slack_token_secret_name})
    slack_token = slack_token_response.payload.data.decode('UTF-8')
    os.environ['SLACK_API_TOKEN'] = slack_token
    print(slack_token)
    return slack_channel


if __name__ == '__main__':
    main()
