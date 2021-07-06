""" Updates the Sample Metadata Database to register new samples """
import csv
import io
import os
import json
import subprocess
from typing import List, Tuple
import click
from google.cloud import storage

from sample_metadata.api.sample_api import SampleApi


def get_csv(bucket_name, prefix) -> Tuple[csv.DictReader, str]:
    """Pull relevant data from csv file in a given bucket
    ==========
    Only one .csv file exists in the upload bucket and previous batch directory.
    """
    client = storage.Client()
    gcs_bucket = client.get_bucket(bucket_name)

    all_blobs = list(client.list_blobs(bucket_name, prefix=prefix))
    csv_path = next(filter(lambda blob: blob.name.endswith('.csv'), all_blobs)).name

    csv_as_text = gcs_bucket.get_blob(csv_path).download_as_text()
    csv_reader = csv.DictReader(io.StringIO(csv_as_text))

    return csv_reader, csv_path


def update_samples(csv_dict_reader, proj):
    """Update the sample status
    When the initial sample list is uploaded (before all sequencing has been uploaded) the status for
    each sample should reflect this i.e. 'New' or 'Registered' or 'Waiting' or something.
    Then, when the sample is uploaded this sample status should change to 'Active'.
    This function will decide what samples should have this status change to upload."""

    # Pull a list containing the sample ID's that don't have gvcfs
    sapi = SampleApi()
    previous_samples_internal = sapi.get_samples_with_gvcfs(proj)  # TODO: Implement
    previous_samples_external = sapi.get_external_ids(
        previous_samples_internal, proj
    )  # TODO: Implement

    # The samples list from the csv file is cumulative.

    # Pull the data from the CSV iterable.

    samples: List[str] = []
    sample_metadata = []

    for sample_dict in csv_dict_reader:
        samples.append(sample_dict['sample.sample_name'])
        sample_metadata.append(sample_dict)

    # Determine the samples in the latest upload.
    latest_upload = samples - previous_samples_external

    sample_meta_map = {d['sample.sample_name']: d for d in sample_metadata}
    for sample in latest_upload:
        metadata = sample_meta_map[sample]
        metadata_json = json.dumps(list(metadata)[0], indent=2)
        sapi.update_metadata(sample, metadata_json)  # TODO: IMPLEMENT.

        sapi.update_status(sample, 'Uploaded')  # TODO: IMPLEMENT.


@click.command()
@click.option('--batch_number')
def upload_samples(batch_number: int):
    """ Update Sample Metadata DB with new samples """

    # Initialise Variables
    project = os.getenv('HAIL_BILLING_PROJECT')
    bucket = f'cpg-{project}-upload'
    metadata_bucket = f'cpg-{project}-main-metadata'
    batch_path = f'batch{batch_number}'

    #  Get the CSV file.
    csv_reader, csv_path = get_csv(bucket, batch_path)

    # Update the status in the DB for the newly uploaded samples.
    update_samples(csv_reader, project)

    # Move the csv when its all done.
    subprocess.run(
        ['gsutil', 'mv', f'gs://{csv_path}', f'gs://{metadata_bucket}/{batch_path}/'],
        check=True,
    )


if __name__ == '__main__':
    upload_samples()  # pylint: disable=E1120
