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
from sample_metadata.apis import SampleApi
from sample_metadata.apis import AnalysisApi
from sample_metadata.apis import SequenceApi
from sample_metadata.models import AnalysisType, AnalysisStatus
from sample_metadata.model.analysis_model import AnalysisModel
from sample_metadata.exceptions import ServiceException
from joint_calling.upload_processor import batch_move_files, SampleGroup


def determine_samples(proj):
    """ Determine which samples should be processed """
    # Pull a list containing the sample ID's that don't have gvcfs and that also
    # Have sequencing metadata attached.
    aapi = AnalysisApi()
    sapi = SampleApi()
    seqapi = SequenceApi()

    samples_without_analysis = aapi.get_all_sample_ids_without_analysis_type(
        'gvcf', proj
    )

    print(samples_without_analysis)

    # print(samples_without_analysis)

    samples_with_sequencing_meta = []

    for sample in samples_without_analysis['sample_ids']:
        print(sample)

        try:
            seq_entry = seqapi.get_sequence_id_from_sample_id(sample, proj)
            print(seq_entry)
            # This will return a dictionary. Check if a dictionary has a key.
            if 'meta' in seq_entry:
                if 'reads' in seq_entry['meta']:
                    samples_with_sequencing_meta.append(sample)

        except ServiceException:
            print(f'Sequencing for {sample} not found')

    print(samples_with_sequencing_meta)

    # Determine the intersection between both lists.
    latest_upload_internal = list(
        set(samples_with_sequencing_meta) & set(samples_without_analysis['sample_ids'])
    )

    latest_upload_external = sapi.get_sample_external_id_map(
        proj, {'internal_ids': latest_upload_internal}
    )

    # latest_upload_internal = sapi.get_sample_id_map(proj, latest_upload_external)

    # latest_upload_external = []

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


def upload_gvcf(sample_group: SampleGroup, proj, main_path):
    """ Creates a new analysis object"""
    sapi = SampleApi()
    aapi = AnalysisApi()
    external_id = {'external_ids': [sample_group.sample_id_external]}
    internal_id_map = sapi.get_sample_id_map(proj, external_id)
    internal_id = list(internal_id_map.values())[0]

    filepath = os.path.join('gs://', main_path, internal_id)

    new_gvcf = AnalysisModel(
        sample_ids=[internal_id],
        type=AnalysisType('gvcf'),
        status=AnalysisStatus('completed'),
        output=filepath,
    )

    aapi.create_new_analysis(proj, new_gvcf)


def upload_cram(sample_group: SampleGroup, proj, archive_path):
    """ Creates a new analysis object"""
    sapi = SampleApi()
    aapi = AnalysisApi()
    external_id = {'external_ids': [sample_group.sample_id_external]}
    internal_id_map = sapi.get_sample_id_map(proj, external_id)
    internal_id = list(internal_id_map.values())[0]

    filepath = os.path.join('gs://', archive_path, internal_id)

    new_gvcf = AnalysisModel(
        sample_ids=[internal_id],
        type=AnalysisType('cram'),
        status=AnalysisStatus('completed'),
        output=filepath,
    )

    aapi.create_new_analysis(proj, new_gvcf)


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
    # metadata_bucket = f'cpg-{project}-main-metadata'
    archive_path = join(f'cpg-{project}-archive', 'cram', batch_path)

    docker_image = os.environ.get('DRIVER_IMAGE')
    key = os.environ.get('GSA_KEY')

    if project is None:
        raise ValueError('HAIL_BILLING_PROJECT must be set')

    # Determine the analysis results (i.e. list of gvcfs and crams) to be moved
    samples_external_ids: List[str] = []  # List of external sample IDs
    main_files: List[SampleGroup] = []
    archive_files: List[SampleGroup] = []

    #  Get the CSV file.
    # csv_reader, csv_path = get_csv(upload_bucket, upload_prefix)

    samples_external_ids = determine_samples(project)

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
            project,
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
        status_job.call(upload_gvcf, sample_group, main_path)
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
            project,
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
        status_job.call(upload_cram, sample_group, archive_path)
        status_job.depends_on(validate_job)
        archive_jobs.append(status_job)

    batch.run()


if __name__ == '__main__':
    # run_processor()  # pylint: disable=E1120
    determine_samples('viviandev')
