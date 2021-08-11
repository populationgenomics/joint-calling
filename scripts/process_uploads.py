#!/usr/bin/env python3
"""
Processes a new batch of samples from the UPLOAD bucket,
into the MAIN and ARCHIVE buckets using the batch_move_files function. 

Assumes that each relevant sample and sample_sequencing has been previously 
uploaded. 
"""
import os
from os.path import join
from typing import List, Optional, Tuple
import hailtop.batch as hb

# Import SampleAPI
from sample_metadata.api import SampleApi
from sample_metadata.api import AnalysisApi
from sample_metadata.api import SequenceApi
from sample_metadata.models.analysis_type import AnalysisType
from sample_metadata.models.analysis_status import AnalysisStatus
from sample_metadata.models.analysis_model import AnalysisModel

# from sample_metadata.exceptions import ServiceException
from joint_calling.upload_processor import batch_move_files, SampleGroup


def determine_samples_archived(proj):
    """ Determine which samples should be processed """
    aapi = AnalysisApi()
    sapi = SampleApi()
    seqapi = SequenceApi()

    samples_without_analysis = aapi.get_all_sample_ids_without_analysis_type(
        'gvcf', proj
    )

    samples_with_sequencing_meta = []

    # Determines which sequences have had their metadata fields updated
    # (This metadata is updated on initial notification of upload)
    for internal_sample_id in samples_without_analysis['sample_ids']:
        seq_id = int(seqapi.get_sequence_id_from_sample_id(internal_sample_id, proj))
        seq_entry = seqapi.get_sequence_by_id(seq_id, proj)
        full_external_id_batch_mapping = []

        # This will return a dictionary. Check if a dictionary has a key.
        if seq_entry['meta'] is not None:
            if 'gvcf' in seq_entry['meta']:
                batch = seq_entry['meta']['gvcf'].get('batch')
                samples_with_sequencing_meta.append(internal_sample_id)
                external_map = sapi.get_sample_id_map_by_internal(
                    proj, [internal_sample_id]
                )
                full_external_id_batch_mapping.append(
                    {
                        'internal_id': internal_sample_id,
                        'external_id': external_map[internal_sample_id],
                        'batch': batch,
                    }
                )

    # Intersection determines the sequencing that is ready to be processed, but has not
    # yet been.
    latest_upload_internal = list(
        set(samples_with_sequencing_meta) & set(samples_without_analysis['sample_ids'])
    )

    latest_upload_external = []

    # Map back to the external IDs (required for move)
    if latest_upload_internal:
        latest_upload_external = sapi.get_sample_id_map_by_internal(
            proj, latest_upload_internal
        )

    latest_upload_external_batch = []

    # Links batch data with external ID
    for sample_id in latest_upload_external:
        mapped_batch = next(
            batch_map
            for batch_map in full_external_id_batch_mapping
            if batch_map['internal_id'] == sample_id
        )

        latest_upload_external_batch.append(mapped_batch)

    return latest_upload_external_batch


def generate_file_list(
    external_sample_ids,
) -> Tuple[List[SampleGroup], List[SampleGroup]]:
    """ Generate list of expected files, given a list of external sample ids"""
    main_files: List[SampleGroup] = []
    archive_files: List[SampleGroup] = []

    for sample_batch_pair in external_sample_ids:
        external_id = sample_batch_pair['external_id']
        internal_id = sample_batch_pair['internal_id']
        batch_id = sample_batch_pair['batch']
        sample_group_main = SampleGroup(
            sample_id_external=external_id,
            sample_id_internal=internal_id,
            data_file=f'{external_id}.g.vcf.gz',
            index_file=f'{external_id}.g.vcf.gz.tbi',
            md5=f'{external_id}.g.vcf.gz.md5',
            batch_number=batch_id,
        )
        main_files.append(sample_group_main)
        sample_group_archive = SampleGroup(
            sample_id_external=external_id,
            sample_id_internal=internal_id,
            data_file=f'{external_id}.cram',
            index_file=f'{external_id}.cram.crai',
            md5=f'{external_id}.cram.md5',
            batch_number=batch_id,
        )
        archive_files.append(sample_group_archive)

    return main_files, archive_files


def determine_samples(proj) -> Tuple[List[SampleGroup], List[SampleGroup]]:
    """ Determine which samples should be processed """
    print('entered here')
    aapi = AnalysisApi()
    sapi = SampleApi()
    seqapi = SequenceApi()

    samples_without_analysis = aapi.get_all_sample_ids_without_analysis_type(
        'gvcf', proj
    )

    sample_ids_without_analysis = samples_without_analysis['sample_ids']
    sequences = seqapi.get_sequences_by_ids(
        sample_ids=sample_ids_without_analysis, project=proj
    )

    external_sample_mapping = sapi.get_sample_id_map_by_internal(
        proj, sample_ids_without_analysis
    )

    main_files = []
    archive_files = []
    print(sequences)

    for seq_entry in sequences:
        print(seq_entry)
        print(type(seq_entry))

        internal_sample_id = seq_entry.get('sample_id')

        if seq_entry['meta'] is not None:
            batch_number = seq_entry['meta'].get('batch')

            if 'gvcf' in seq_entry['meta']:
                data_file = seq_entry['meta']['gvcf'].get('location')
                # Note this will handle only one secondary file.
                index_file = seq_entry['meta']['gvcf']['secondaryFiles'][0].get(
                    'location'
                )

                sample_group_main = SampleGroup(
                    sample_id_external=external_sample_mapping.get(internal_sample_id),
                    sample_id_internal=internal_sample_id,
                    data_file=data_file,
                    index_file=index_file,
                    md5=data_file + '.md5',
                    batch_number=batch_number,
                )

                main_files.append(sample_group_main)

            if 'reads' in seq_entry['meta']:

                for read in seq_entry['meta']['reads']:
                    data_file = read.get('location')
                    index_file = read['secondaryFiles'][0].get('location')
                    sample_group_archive = SampleGroup(
                        sample_id_external=external_sample_mapping.get(
                            internal_sample_id
                        ),
                        sample_id_internal=internal_sample_id,
                        data_file=data_file,
                        index_file=index_file,
                        md5=data_file + '.md5',
                        batch_number=batch_number,
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

    if not python_job:
        job.command(
            'gcloud -q auth activate-service-account --key-file=/gsa-key/key.json'
        )

    return job


def create_analysis_in_sm_db(sample_group: SampleGroup, proj, path, analysis_type):
    """ Creates a new analysis object"""
    aapi = AnalysisApi()

    internal_id = sample_group.sample_id_internal

    if analysis_type == 'gvcf':
        file_extension = '.g.vcf.gz'
    else:
        file_extension = '.cram'

    # TODO Rebuild this file path
    filepath = os.path.join('gs://', path, internal_id + file_extension)

    new_analysis = AnalysisModel(
        sample_ids=[internal_id],
        type=analysis_type,
        status=AnalysisStatus.COMPLETED,
        output=filepath,
    )

    aapi.create_new_analysis(proj, new_analysis)


def validate_md5(
    job: hb.batch.job, sample_group: SampleGroup, upload_path: str
) -> hb.batch.job:
    """Each gvcf and cram file is uploaded with a corresponding md5 checksum.
    This function appends commands to the job that moves these files, to validate each .md5 file."""

    # TODO: Re-work this function.

    file_extension = sample_group.data_file[len(sample_group.sample_id_internal) :]
    file_path = f'{sample_group.sample_id_internal}{file_extension}'

    # Generate paths to files that are being validated
    path_to_data = join(
        'gs://',
        upload_path,
        f'batch{sample_group.batch_number}',
        file_path,
    )
    # extension = sample_group.data_file[len(sample)]

    path_to_md5 = join(
        'gs://',
        upload_path,
        f'batch{sample_group.batch_number}',
        f'{file_path}.md5',
    )

    # Calculate md5 checksum.
    job.command(
        f'gsutil cat {path_to_data} | md5sum | sed "s/-/{sample_group.data_file}/" > /tmp/{sample_group.md5}'
    )

    job.command(
        f'diff <(cat /tmp/{sample_group.md5} | cut -d " " -f1 ) <(gsutil cat {path_to_md5} | cut -d " " -f1 )'
    )

    return job


def run_processor():
    """Set up and execute batch_move_files, validate_md5 and update_status"""

    # Setting up inputs for batch_move_files
    project = os.getenv('HAIL_BILLING_PROJECT')
    upload_path = join(f'cpg-{project}-main-upload')
    main_bucket = f'cpg-{project}-main'
    main_path = join(main_bucket, 'gvcf')
    archive_path = join(f'cpg-{project}-archive', 'cram')

    sm_project = 'viviandev'
    # sm_project = project.replace('-', '')
    docker_image = os.environ.get('DRIVER_IMAGE')
    key = os.environ.get('GSA_KEY')

    if project is None:
        raise ValueError('HAIL_BILLING_PROJECT must be set')

    # Determine the analysis results (i.e. list of gvcfs and crams) to be moved
    # samples_external_ids: List[str] = []  # List of external sample IDs
    main_files: List[SampleGroup] = []
    archive_files: List[SampleGroup] = []

    # samples_external_ids = determine_samples(sm_project)

    # main_files, archive_files = generate_file_list(samples_external_ids)

    main_files, archive_files = determine_samples(sm_project)
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
        full_path = os.path.join(main_path, f'batch{sample_group.batch_number}')
        status_job.call(
            create_analysis_in_sm_db,
            sample_group,
            sm_project,
            full_path,
            AnalysisType.GVCF,
        )
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
        status_job.call(
            create_analysis_in_sm_db,
            sample_group,
            sm_project,
            archive_path,
            AnalysisType.CRAM,
        )
        status_job.depends_on(validate_job)
        archive_jobs.append(status_job)

    batch.run()


if __name__ == '__main__':
    run_processor()  # pylint: disable=E1120
