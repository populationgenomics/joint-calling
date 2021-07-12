#!/usr/bin/env python3
"""
Processes a new batch of samples from the UPLOAD bucket,
into the MAIN and ARCHIVE buckets using the batch_move_files function. 

Assumes that the csv file has been processed in a previous step. 
"""
import csv
import io
import os
import json
from os.path import join
from typing import List, Optional, Tuple
import click
import hailtop.batch as hb
from google.cloud import storage

# Import SampleAPI
from sample_metadata.api.sample_api import SampleApi
from joint_calling.upload_processor import batch_move_files, SampleGroup


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
    previous_samples_internal = sapi.get_all_sample_ids_without_analysis_type(
        proj, 'gvcf'
    )
    previous_samples_external = sapi.get_external_ids(
        previous_samples_internal, proj
    )  # TODO: Implement

    samples: List[str] = []
    sample_metadata = []

    for sample_dict in csv_dict_reader:
        samples.append(sample_dict['sample.sample_name'])
        sample_metadata.append(sample_dict)

    # Determine the samples in the latest upload.
    latest_upload_external = samples - previous_samples_external
    latest_upload_internal = sapi.get_internal_ids(latest_upload_external, proj)

    sample_meta_map = {d['sample.sample_name']: d for d in sample_metadata}
    for sample in latest_upload_internal:
        metadata = sample_meta_map[sample]
        metadata_json = json.dumps(list(metadata)[0], indent=2)

        sapi.create_new_analysis(sample, 'gvcf')  # TODO, Fix inputs
        sapi.create_new_analysist(sample, 'cram')

        sapi.update_metadata(sample, metadata_json)  # TODO: IMPLEMENT.

    return latest_upload_external


def generate_file_list(
    external_sample_ids: List[str],
) -> Tuple[List[SampleGroup], List[SampleGroup]]:
    """ Generate list of expected files, given a list of external sample ids"""
    main_files: List[SampleGroup] = []
    archive_files: List[SampleGroup] = []

    for external_id in external_sample_ids:
        sample_group_main = SampleGroup(
            external_id,
            f'{external_id}.g.vcf.gz',
            f'{external_id}.g.vcf.gz.tbi',
            f'{external_id}.g.vcf.gz.md5',
        )
        main_files.append(sample_group_main)
        sample_group_archive = SampleGroup(
            external_id,
            f'{external_id}.cram',
            f'{external_id}.cram.crai',
            f'{external_id}.cram.md5',
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


def update_status(sample_group: SampleGroup):
    """ Updates the status of a SampleSequence """
    sapi = SampleApi()
    sapi.update_sequencing_status_from_external_sample_id(
        sample_group.sample_id_external, 'Upload Successful'
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
    metadata_bucket = f'cpg-{project}-main-metadata'
    archive_path = join(f'cpg-{project}-archive', 'cram', batch_path)

    docker_image = os.environ.get('DRIVER_IMAGE')
    key = os.environ.get('GSA_KEY')

    # Determine the analysis results (i.e. list of gvcfs and crams) to be moved
    samples_external_ids: List[str] = []  # List of external sample IDs
    main_files: List[SampleGroup] = []
    archive_files: List[SampleGroup] = []

    #  Get the CSV file.
    csv_reader, csv_path = get_csv(upload_bucket, upload_prefix)

    samples_external_ids = create_analysis(csv_reader, project)

    # samples = get_samples_from_db(project)
    main_files, archive_files = generate_file_list(samples_external_ids)

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
        main_jobs.append(status_job)

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
        archive_jobs.append(status_job)

    # Move the csv file to the metadata bucket, after the gVCFs and CRAMs have been
    # moved to the main and archive buckets.
    csv_job = setup_job(batch, f'Move {csv_path}', docker_image)
    csv_job.command(f'gsutil mv gs://{csv_path} gs://{metadata_bucket}/{batch_path}/')
    csv_job.depends_on(*main_jobs)
    csv_job.depends_on(*archive_jobs)

    batch.run()


if __name__ == '__main__':
    run_processor()  # pylint: disable=E1120
