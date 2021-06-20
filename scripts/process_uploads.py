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
from typing import List, Optional
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
    upload_bucket: str, upload_prefix: str, metadata_bucket: str, metadata_prefix: str
):
    """Determine files that should be moved to main vs archive bucket.
    Determines the difference between the latest CSV file within upload
    and the CSV file from the most recent batch as the list of samples to
    be moved.
    """
    all_samples: List[str] = []
    previous_samples: List[str] = []

    all_samples, csv_path_upload = samples_from_csv(upload_bucket, upload_prefix)
    previous_samples, _ = samples_from_csv(metadata_bucket, metadata_prefix)

    samples = set(all_samples) - set(previous_samples)

    curr_csv_file_path = join(upload_bucket, csv_path_upload)

    return samples, curr_csv_file_path


def generate_file_list(samples: List[str]):
    """ Generate list of expected files, given a list of sampleIDs """
    main_files = []
    archive_files = []
    for s in samples:
        sample_group_main = [f'{s}.g.vcf.gz', f'{s}.g.vcf.gz.tbi', f'{s}.g.vcf.gz.md5']
        main_files.append(sample_group_main)
        sample_group_archive = [f'{s}.cram', f'{s}.cram.crai', f'{s}.cram.md5']
        archive_files.append(sample_group_archive)

    return main_files, archive_files


def setup_job(batch: hb.batch, name: str, docker_image: Optional[str]):
    """ Define a new hail batch job """

    new_j = batch.new_job(name=name)
    if docker_image is not None:
        new_j.image(docker_image)
    # Uncomment for local testing.
    # key = os.environ.get('GSA_KEY')
    # new_j.command(f"echo '{key}' > /tmp/key.json")
    # new_j.command(f'gcloud -q auth activate-service-account --key-file=/tmp/key.json')
    new_j.command(
        'gcloud -q auth activate-service-account --key-file=/gsa-key/key.json'
    )

    return new_j


def validate_md5(job: hb.batch.job, sample_group: List[str], upload_path: str):
    """Each gvcf and cram file is uploaded with a corresponding md5 checksum.
    This function validates each .md5 file, by recalculating the md5 checksum
    and comparing this with the uploaded .md5."""

    # Generate paths to files that are being validated
    sample_gvcf = sample_group[0]
    md5_file = sample_group[2]
    path_to_gvcf = join('gs://', upload_path, sample_gvcf)
    path_to_md5 = join('gs://', upload_path, md5_file)

    # Calculate md5 checksum.
    job.command(
        f'gsutil cat {path_to_gvcf} | md5sum | sed "s/-/{sample_gvcf}/" > {job.ofile}'
    )
    job.command(f'diff <( cat {job.ofile} ) <(gsutil cat {path_to_md5})')

    return job


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
    metadata_bucket = f'cpg-{project}-main-metadata'
    archive_path = join(f'cpg-{project}-archive', 'cram', batch_path)

    docker_image = os.environ.get('DRIVER_IMAGE')
    key = os.environ.get('GSA_KEY')

    # Determine files to be processed
    samples, csv_path = determine_samples(
        upload_bucket, upload_prefix, metadata_bucket, prev_batch
    )
    main_files, archive_files = generate_file_list(samples)

    service_backend = hb.ServiceBackend(
        billing_project=project,
        bucket=os.getenv('HAIL_BUCKET'),
    )

    batch = hb.Batch(name='Process files', backend=service_backend)

    main_jobs = []
    for sample_group in main_files:

        # Moving the files to the main bucket
        sample_group_main_jobs = batch_move_files(
            batch,
            sample_group,
            upload_path,
            main_path,
            docker_image,
            key,
        )

        validate_job = setup_job(batch, 'Validate MD5', docker_image)
        validate_job = validate_md5(validate_job, sample_group, main_path)
        validate_job.depends_on(*sample_group_main_jobs)
        main_jobs.append(validate_job)

    archive_jobs = []
    for sample_group in archive_files:
        # Moving the files to the archive bucket
        sample_group_archive_jobs = batch_move_files(
            batch,
            sample_group,
            upload_path,
            archive_path,
            docker_image,
            key,
        )
        validate_job = setup_job(batch, 'Validate MD5', docker_image)
        validate_job = validate_md5(validate_job, sample_group, main_path)
        validate_job.depends_on(*sample_group_archive_jobs)
        archive_jobs.append(validate_job)

    # Move the csv file to the metadata bucket, after the gVCFs and CRAMs have been
    # moved to the main and archive buckets.
    csv_job = setup_job(batch, f'Move {basename(csv_path)}', docker_image)
    csv_job.command(f'gsutil mv gs://{csv_path} gs://{metadata_bucket}/{batch_path}/')
    csv_job.depends_on(*main_jobs)
    csv_job.depends_on(*archive_jobs)

    batch.run()


if __name__ == '__main__':
    run_processor()  # pylint: disable=E1120
