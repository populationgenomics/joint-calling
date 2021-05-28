#!/usr/bin/env python3
"""
Processes files uploaded to gcp
"""

import os
import csv
from os.path import join, basename
import subprocess
import tempfile
import shutil
from typing import List
import hail as hl
import hailtop.batch as hb
from joint_calling.upload_processor import batch_move_files


def samples_from_csv(local_csv_path):
    """ Determines list of samples from a given csv file. """

    samples: List[str] = []
    with open(local_csv_path) as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            samples.append(row['sample.sample_name'])

    return samples


def determine_samples(upload_bucket: str, previous_batch_path: str):
    """Determine files that should be moved to main vs archive bucket
    Parameters
    =========
    upload_bucket : Name of the upload bucket
    previous_batch_path : Previous batch folder

    Assumptions
    ==========
    Only one .csv file exists in the bucket.
    """
    local_tmp_dir = tempfile.mkdtemp()
    all_samples: List[str] = []
    previous_samples: List[str] = []

    # Pull the samples listed in the most recent CSV file
    cmd = f'gsutil ls \'gs://{upload_bucket}/*.csv\''
    curr_csv_file_path = subprocess.check_output(cmd, shell=True).decode().strip()
    local_curr_csv_path = join(local_tmp_dir, basename(curr_csv_file_path))
    subprocess.run(
        f'gsutil cp {curr_csv_file_path} {local_curr_csv_path}', check=False, shell=True
    )

    all_samples = samples_from_csv(local_curr_csv_path)

    # Pull the samples listed in the CSV file from the previous batch
    cmd = f'gsutil ls \'gs://{previous_batch_path}/*.csv\''
    prev_csv_file_path = subprocess.check_output(cmd, shell=True).decode().strip()
    local_prev_csv_path = join(local_tmp_dir, basename(prev_csv_file_path))
    subprocess.run(
        f'gsutil cp {prev_csv_file_path} {local_prev_csv_path}', check=False, shell=True
    )

    previous_samples = samples_from_csv(local_prev_csv_path)

    samples = set(all_samples) - set(previous_samples)

    shutil.rmtree(local_tmp_dir)

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
    e.g. "batch0"

    """

    # Setting up inputs for batch_move_files
    project = os.getenv('HAIL_BILLING_PROJECT')
    upload_prefix = os.path.join(f'cpg-{project}-upload')
    main_prefix = os.path.join(f'cpg-{project}-main', 'gvcf', batch_number)
    prev_prefix = os.path.join(f'cpg-{project}-main', 'gvcf', prev_batch)
    docker_image = os.environ.get('DRIVER_IMAGE')
    key = os.environ.get('GSA_KEY')

    # Determine files to be processed
    samples, csv_path = determine_samples(upload_prefix, prev_prefix)
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

    archive_prefix = os.path.join(f'cpg-{project}-archive', 'cram', batch_number)

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
