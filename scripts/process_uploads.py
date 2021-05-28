#!/usr/bin/env python3
"""
Processes files uploaded to gcp
"""

import os
import io
import csv
from os.path import join, basename
import subprocess
from typing import List
import hail as hl
import hailtop.batch as hb
from google.cloud import storage
from joint_calling.upload_processor import batch_move_files


def samples_from_csv(bucket_name, path):
    """ Determines list of samples from a given csv file. """

    client = storage.Client()
    bucket = client.get_bucket(bucket_name)
    # full_path = os.path.join('gs://', bucket, path)

    cmd = f'gsutil ls \'gs://{bucket_name}/{path}\''
    csv_path = subprocess.check_output(cmd, shell=True).decode().strip()
    blob = bucket.get_blob(basename(csv_path))
    data = blob.download_as_string().decode('utf-8')
    csv_reader = csv.DictReader(io.StringIO(data))

    samples: List[str] = []
    for row in csv_reader:
        samples.append(row['sample.sample_name'])

    return samples


def determine_samples(upload_bucket: str, main_bucket: str, previous_batch_path: str):
    """Determine files that should be moved to main vs archive bucket.
    Determines the difference between the latest CSV file within upload
    and the CSV file from the most recent batch as the list of samples to
    be moved.

    Assumptions
    ==========
    Only one .csv file exists in the upload bucket and previous batch directory.
    """
    # local_tmp_dir = tempfile.mkdtemp()
    all_samples: List[str] = []
    previous_samples: List[str] = []

    # cmd = f'gsutil ls \'gs://{upload_bucket}/*.csv\''
    # curr_csv_file_path = subprocess.check_output(cmd, shell=True).decode().strip()
    # blob = bucket.get_blob(basename(curr_csv_file_path))
    # data = blob.download_as_string().decode("utf-8")
    # reader_list = csv.DictReader(io.StringIO(data))
    # upload_path = os.path.join('gs://', upload_bucket, upload_path)

    # Pull the samples listed in the most recent CSV file
    cmd = f'gsutil ls \'gs://{upload_bucket}\''
    curr_csv_file_path = subprocess.check_output(cmd, shell=True).decode().strip()
    # local_curr_csv_path = join(local_tmp_dir, basename(curr_csv_file_path))
    # subprocess.run(
    #     f'gsutil cp {curr_csv_file_path} {local_curr_csv_path}', check=False, shell=True
    # )

    all_samples = samples_from_csv(upload_bucket, '')

    # Pull the samples listed in the CSV file from the previous batch
    # cmd = f'gsutil ls \'gs://{previous_batch_path}/*.csv\''
    # prev_csv_file_path = subprocess.check_output(cmd, shell=True).decode().strip()
    # blob = bucket.get_blob(basename(curr_csv_file_path))
    # data = blob.download_as_string().decode("utf-8")
    # reader_list = csv.DictReader(io.StringIO(data))

    # local_prev_csv_path = join(local_tmp_dir, basename(prev_csv_file_path))
    # subprocess.run(
    #     f'gsutil cp {prev_csv_file_path} {local_prev_csv_path}', check=False, shell=True
    # )

    previous_samples = samples_from_csv(main_bucket, previous_batch_path)

    samples = set(all_samples) - set(previous_samples)

    # shutil.rmtree(local_tmp_dir)

    return samples, curr_csv_file_path


def generate_file_list(samples: List[str]):
    """ Generate list of expected files, given a list of sampleIDs """
    main_files: List[str] = []
    archive_files: List[str] = []
    for s in samples:
        main_files.extend([f'{s}.g.vcf.gz', f'{s}.g.vcf.gz.tbi', f'{s}.g.vcf.gz.md5'])
        archive_files.extend([f'{s}.cram', f'{s}.cram.crai', f'{s}.cram.md5'])

    return main_files, archive_files


def run_processor(batch_number: str, prev_batch: str):
    """Set up and execute batch_move_files

    Parameters
    =========
    batch_number: str
    Corresponds to the folder within main that the files should be moved to
    e.g. "batch2"

    prev_batch: str
    Corresponds to the most recent batch subdirectory.
    e.g. "batch1"

    """

    # Setting up inputs for batch_move_files
    project = os.getenv('HAIL_BILLING_PROJECT')
    # upload_bucket = os.path.join(f'cpg-{project}-upload')
    # upload_prefix = os.path.join(upload_bucket)
    # main_bucket = os.path.join(f'cpg-{project}-main')
    # main_prefix = os.path.join(main_bucket, 'gvcf', batch_number)
    # prev_prefix = os.path.join('gvcf', prev_batch)
    # archive_prefix = os.path.join(f'cpg-{project}-archive', 'cram', batch_number)

    upload_bucket = f'cpg-{project}-temporary'
    upload_prefix = os.path.join(upload_bucket, 'vivian-test', 'upload')
    main_bucket = f'cpg-{project}-temporary'
    main_prefix = os.path.join(main_bucket, 'vivian-test', 'main', 'gvcf', batch_number)
    prev_prefix = os.path.join('gvcf', prev_batch)
    archive_prefix = os.path.join(
        f'cpg-{project}-temporary', 'vivian-test', 'archive', 'cram', batch_number
    )

    docker_image = os.environ.get('DRIVER_IMAGE')
    key = os.environ.get('GSA_KEY')

    # Determine files to be processed
    samples, csv_path = determine_samples(upload_bucket, main_bucket, prev_prefix)
    main_files, archive_files = generate_file_list(samples)

    # Initialize the service backend.
    hl.init()

    service_backend = hb.ServiceBackend(
        billing_project=project,
        bucket=os.getenv('HAIL_BUCKET'),
    )

    batch = hb.Batch(name='Process files', backend=service_backend)

    # Moving the files to the main bucket
    main_jobs = batch_move_files(
        batch,
        main_files,
        upload_prefix,
        main_prefix,
        docker_image,
        key,
    )

    # Moving the files to the archive bucket
    archive_jobs = batch_move_files(
        batch,
        archive_files,
        upload_prefix,
        archive_prefix,
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
    final_csv_location = join('gs://', main_prefix, basename(csv_path))
    csv_job.command(f'gsutil mv {csv_path} {final_csv_location}')

    batch.run()


if __name__ == '__main__':

    # Run the processor on the 4th batch, i.e. batch3.
    run_processor('batch3', 'batch2')
