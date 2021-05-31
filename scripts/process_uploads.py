#!/usr/bin/env python3
"""
Processes a new batch of samples from the UPLOAD bucket,
into the MAIN and ARCHIVE buckets using the batch_move_files function. 

Determines the files to be moved through a comparison of the 
most recently updated CSV and the CSV uploaded in a previous batch run. 
"""

import os
import io
import csv
from os.path import join, basename
from typing import List
import click
import hailtop.batch as hb
from google.cloud import storage
from joint_calling.upload_processor import batch_move_files


def samples_from_csv(bucket_name, prefix):
    """Determines list of samples from a csv file within a specified bucket

    Assumptions
    ==========
    Only one .csv file exists in the upload bucket and previous batch directory.
    """

    client = storage.Client()
    bucket = client.get_bucket(bucket_name)

    all_blobs = list(client.list_blobs(bucket_name, prefix=prefix))
    csv_path = next(filter(lambda blob: blob.name.endswith('.csv'), all_blobs)).name

    blob = bucket.get_blob(csv_path).download_as_text()
    csv_reader = csv.DictReader(io.StringIO(blob))

    samples: List[str] = []
    for row in csv_reader:
        samples.append(row['sample.sample_name'])

    return samples, csv_path


def determine_samples(
    upload_bucket: str, upload_prefix: str, main_bucket: str, previous_batch_prefix: str
):
    """Determine files that should be moved to main vs archive bucket.
    Determines the difference between the latest CSV file within upload
    and the CSV file from the most recent batch as the list of samples to
    be moved.
    """
    all_samples: List[str] = []
    previous_samples: List[str] = []

    all_samples, csv_path_upload = samples_from_csv(upload_bucket, upload_prefix)
    previous_samples, _ = samples_from_csv(main_bucket, previous_batch_prefix)

    samples = set(all_samples) - set(previous_samples)

    curr_csv_file_path = join(upload_bucket, csv_path_upload)

    return samples, curr_csv_file_path


def generate_file_list(samples: List[str]):
    """ Generate list of expected files, given a list of sampleIDs """
    main_files: List[str] = []
    archive_files: List[str] = []
    for s in samples:
        main_files.extend([f'{s}.g.vcf.gz', f'{s}.g.vcf.gz.tbi', f'{s}.g.vcf.gz.md5'])
        archive_files.extend([f'{s}.cram', f'{s}.cram.crai', f'{s}.cram.md5'])

    return main_files, archive_files


@click.command()
@click.option('--batch_number')
def run_processor(
    batch_number: int,
):
    """Set up and execute batch_move_files

    Parameters
    =========
    batch_number: int
    Indexed starting at 1. Used to navigate subdirectory structure
    e.g. 1 --> batch1

    """

    # Setting up inputs for batch_move_files
    project = os.getenv('HAIL_BILLING_PROJECT')
    batch_path = f'batch{batch_number}'
    prev_batch = f'batch{int(batch_number)-1}'
    upload_bucket = f'cpg-{project}-upload'
    upload_prefix = ''
    upload_path = join(upload_bucket, upload_prefix)
    main_bucket = f'cpg-{project}-main'
    main_prefix = join('gvcf', batch_path)
    main_path = join(main_bucket, main_prefix)
    prev_prefix = join('gvcf', prev_batch)
    archive_path = join(f'cpg-{project}-archive', 'cram', batch_path)

    docker_image = os.environ.get('DRIVER_IMAGE')
    key = os.environ.get('GSA_KEY')

    # Determine files to be processed
    samples, csv_path = determine_samples(
        upload_bucket, upload_prefix, main_bucket, prev_prefix
    )
    main_files, archive_files = generate_file_list(samples)

    service_backend = hb.ServiceBackend(
        billing_project=project,
        bucket=os.getenv('HAIL_BUCKET'),
    )

    batch = hb.Batch(name='Process files', backend=service_backend)

    # Moving the files to the main bucket
    main_jobs = batch_move_files(
        batch,
        main_files,
        upload_path,
        main_path,
        docker_image,
        key,
    )

    # Moving the files to the archive bucket
    archive_jobs = batch_move_files(
        batch,
        archive_files,
        upload_path,
        archive_path,
        docker_image,
        key,
    )

    # Move the csv file, after main and archive buckets have been processed.
    csv_job = batch.new_job(name=f'Move {basename(csv_path)}')
    csv_job.image(docker_image)
    csv_job.command(
        'gcloud -q auth activate-service-account --key-file=/gsa-key/key.json'
    )
    all_jobs = main_jobs + archive_jobs
    csv_job.depends_on(*all_jobs)
    current_csv_location = join('gs://', csv_path)
    final_csv_location = join('gs://', main_path, basename(csv_path))
    csv_job.command(f'gsutil mv {current_csv_location} {final_csv_location}')

    batch.run()


if __name__ == '__main__':
    run_processor()  # pylint: disable=E1120
