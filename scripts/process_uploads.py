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


def filter_files_manual():
    """Determine files that should be moved to main vs archive bucket
    when given a list of file paths"""
    files_for_main = []
    files_for_archive = []
    prefix = 'gs://cpg-tob-wgs-upload/'
    with open('test/data/tob_wgs_batch2.txt') as f:
        for file_path in f:
            file_path = file_path.strip()
            if file_path.endswith(
                ('.g.vcf.gz', '.g.vcf.gz.tbi', '.g.vcf.gz.md5', 'csv')
            ):
                files_for_main.append(file_path[len(prefix) :])
            else:
                files_for_archive.append(file_path[len(prefix) :])

    return files_for_main, files_for_archive


def determine_samples(upload_bucket: str):
    """Determine files that should be moved to main vs archive bucket
    Parameters
    =========
    input_bucket:  List[str] The bucket to be processed

    Assumptions
    ==========
    Only one .csv file exists in the bucket.
    """
    local_tmp_dir = tempfile.mkdtemp()
    samples: List[str] = []

    # Pull the CSV file that contains the list of samples in the current batch

    cmd = f'gsutil ls \'gs://{upload_bucket}/*.csv\''
    csv_file_path = subprocess.check_output(cmd, shell=True).decode().strip()
    local_csv_path = join(local_tmp_dir, basename(csv_file_path))
    subprocess.run(
        f'gsutil cp {csv_file_path} {local_csv_path}', check=False, shell=True
    )

    with open(local_csv_path) as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            samples.append(row['sample.sample_name'])

    shutil.rmtree(local_tmp_dir)
    return samples


def generate_file_list(samples: List[str]):
    """ Generate list of expected files, given a list of sampleIDs """
    main_files: List[str] = []
    archive_files: List[str] = []
    for s in samples:
        main_files.extend([f'{s}.g.vcf.gz', f'{s}.g.vcf.gz.tbi', f'{s}.g.vcf.gz.md5'])
        archive_files.extend([f'{s}.cram', f'{s}.cram.crai', f'{s}.cram.md5'])

    return main_files, archive_files


def run_processor():
    """ Execute upload processor """

    # Setting up inputs for batch_move_files
    project = os.getenv('HAIL_BILLING_PROJECT')
    upload_prefix = os.path.join(f'cpg-{project}-upload')
    main_prefix = os.path.join(f'cpg-{project}', 'gvcf', 'batch2')
    docker_image = os.environ.get('DRIVER_IMAGE')
    key = os.environ.get('GSA_KEY')

    # Process input file of sample names
    samples = determine_samples(upload_prefix)
    main_files, archive_files = generate_file_list(samples)

    # Initialize the service backend.
    hl.init()

    service_backend = hb.ServiceBackend(
        billing_project=project,
        bucket=os.getenv('HAIL_BUCKET'),
    )

    batch = hb.Batch(name='Process files', backend=service_backend)

    # Moving the files to the main bucket
    batch_move_files(
        batch,
        main_files,
        upload_prefix,
        main_prefix,
        docker_image,
        key,
    )

    archive_prefix = os.path.join(f'cpg-{project}-archive', 'cram', 'batch2')

    # Moving the files to the archive bucket
    batch_move_files(
        batch,
        archive_files,
        upload_prefix,
        archive_prefix,
        docker_image,
        key,
    )

    batch.run()
