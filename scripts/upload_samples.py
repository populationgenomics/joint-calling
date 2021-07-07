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


def create_analysis(csv_dict_reader, proj):
    """New analysis objects created for the gvcf and cram for each sample
    Assumptions:
    SampleSequence object created previously.
    We are exclusively processing crams and gvcfs.
    Samples with gvcfs have crams, and samples without gvcfs don't.
    """

    # Pull a list containing the sample ID's that don't have gvcfs
    sapi = SampleApi()
    previous_samples_internal = sapi.get_samples_with_gvcfs(proj)  # TODO: Implement
    previous_samples_external = sapi.get_external_ids(
        previous_samples_internal, proj
    )  # TODO: Implement

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

        sapi.create_analysis_object(sample, 'gvcf')
        sapi.create_analysis_object(sample, 'cram')

        sapi.update_metadata(sample, metadata_json)  # TODO: IMPLEMENT.

        sapi.update_sequencing_status(sample, 'Ready')  # TODO: IMPLEMENT.


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

    # Create new analysis objects for each sequencing result
    create_analysis(csv_reader, project)

    # Move the csv when its all done.
    subprocess.run(
        ['gsutil', 'mv', f'gs://{csv_path}', f'gs://{metadata_bucket}/{batch_path}/'],
        check=True,
    )


if __name__ == '__main__':
    upload_samples()  # pylint: disable=E1120
