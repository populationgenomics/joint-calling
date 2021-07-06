#!/usr/bin/env python3
"""
Processes a new batch of samples from the UPLOAD bucket,
into the MAIN and ARCHIVE buckets using the batch_move_files function. 

Assumes that the csv file has been processed in a previous step. 
"""

import os
from os.path import join
from typing import List, Optional, Tuple
import click
import hailtop.batch as hb

# Import SampleAPI
from sample_metadata.api.sample_api import SampleApi
from joint_calling.upload_processor import batch_move_files, SampleGroup


# Reading CSV pulled out into separate function.


def get_samples_from_db(proj) -> List[str]:
    """ Pulls list of samples to be moved from SampleMetadata DB """

    sapi = SampleApi()
    samples: List[str] = []

    samples = sapi.get_analysis_with_status(proj, 'gvcf', 'Ready')  # TODO: Implement
    samples_external_ids = sapi.get_external_ids(samples, proj)

    return samples_external_ids


def generate_file_list(
    samples: List[str],
) -> Tuple[List[SampleGroup], List[SampleGroup]]:
    """ Generate list of expected files, given a list of samples """
    main_files = []
    archive_files = []
    for s in samples:
        sample_group_main = SampleGroup(
            s, f'{s}.g.vcf.gz', f'{s}.g.vcf.gz.tbi', f'{s}.g.vcf.gz.md5'
        )
        main_files.append(sample_group_main)
        sample_group_archive = SampleGroup(
            s, f'{s}.cram', f'{s}.cram.crai', f'{s}.cram.md5'
        )
        archive_files.append(sample_group_archive)

    return main_files, archive_files


def setup_job(
    batch: hb.batch, name: str, docker_image: Optional[str], python_job: bool = False
) -> hb.batch.job:
    """ Returns a new Hail Batch job that activates the Google service account. """

    job = batch.new_python_job(name=name) if python_job else batch.new_job(name=name)

    if docker_image is not None:
        job.image(docker_image)

    job.command('gcloud -q auth activate-service-account --key-file=/gsa-key/key.json')

    return job


# Optional.
def update_status(sample_group: SampleGroup):
    """ Updates the status of a given sample """
    # Update the status of the sample group that you just uploaded.
    # Condition here depending on the file type.
    sapi = SampleApi()

    sapi.update_sequencing_status(
        sample_group.sample_id_external, 'Uploaded'
    )  # TO IMPLEMENT: API CALL


def validate_md5(
    job: hb.batch.job, sample_group: SampleGroup, upload_path: str
) -> hb.batch.job:
    """Each gvcf and cram file is uploaded with a corresponding md5 checksum.
    This function appends commands to the job that moves these files, to validate each .md5 file."""

    # Generate paths to files that are being validated
    path_to_data = join('gs://', upload_path, sample_group.data_file)
    path_to_md5 = join('gs://', upload_path, sample_group.md5)

    # Calculate md5 checksum.
    job.command(
        f'gsutil cat {path_to_data} | md5sum | sed "s/-/{sample_group.data_file}/" > /tmp/{sample_group.md5}'
    )
    job.command(f'diff /tmp/{sample_group.md5} <(gsutil cat {path_to_md5})')

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
    upload_bucket = f'cpg-{project}-upload'
    upload_prefix = ''
    upload_path = join(upload_bucket, upload_prefix)
    main_bucket = f'cpg-{project}-main'
    main_prefix = join('gvcf', batch_path)
    main_path = join(main_bucket, main_prefix)
    archive_path = join(f'cpg-{project}-archive', 'cram', batch_path)

    docker_image = os.environ.get('DRIVER_IMAGE')
    key = os.environ.get('GSA_KEY')

    samples = get_samples_from_db(project)
    main_files, archive_files = generate_file_list(samples)

    service_backend = hb.ServiceBackend(
        billing_project=project,
        bucket=os.getenv('HAIL_BUCKET'),
    )

    batch = hb.Batch(name='Process files', backend=service_backend)

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

        validate_job = setup_job(
            batch, f'Validate MD5 {sample_group.sample_id_external}', docker_image
        )
        validate_job = validate_md5(validate_job, sample_group, main_path)
        validate_job.depends_on(*sample_group_main_jobs)

        # After each sample is moved, update the status in the SampleMetadata DB to reflect the change in status
        status_job = setup_job(
            batch,
            f'Update Status {sample_group.sample_id_external}',
            docker_image,
            True,
        )
        status_job.call(update_status, sample_group)
        status_job.depends_on(validate_job)

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
        validate_job = setup_job(
            batch, f'Validate MD5 {sample_group.sample_id_external}', docker_image
        )
        validate_job = validate_md5(validate_job, sample_group, archive_path)
        validate_job.depends_on(*sample_group_archive_jobs)

        status_job = setup_job(
            batch,
            f'Update Status {sample_group.sample_id_external}',
            docker_image,
            True,
        )
        status_job.call(update_status, sample_group)
        status_job.depends_on(validate_job)

    batch.run()


if __name__ == '__main__':
    run_processor()  # pylint: disable=E1120
