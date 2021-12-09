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

from sample_metadata.apis import SampleApi, AnalysisApi, SequenceApi
from sample_metadata.models import AnalysisType, AnalysisStatus, AnalysisModel

from joint_calling.upload_processor import batch_move_files, SampleGroup, FileGroup


def determine_samples(proj) -> Tuple[List[SampleGroup], List[SampleGroup]]:
    """ Determine which samples should be processed """
    aapi = AnalysisApi()
    sapi = SampleApi()
    seqapi = SequenceApi()

    samples_without_analysis = aapi.get_all_sample_ids_without_analysis_type(
        'gvcf', proj
    )

    # Take the first 300 samples without analysis
    sample_ids_without_analysis = samples_without_analysis['sample_ids'][:300]
    sequences = seqapi.get_sequences_by_sample_ids(
        request_body=sample_ids_without_analysis
    )

    external_sample_mapping = sapi.get_sample_id_map_by_internal(
        request_body=sample_ids_without_analysis
    )

    main_files = []
    archive_files = []

    for seq_entry in sequences:
        internal_sample_id = seq_entry.get('sample_id')

        if seq_entry['meta'] is not None:
            batch_number = seq_entry['meta'].get('batch')

            if 'gvcf' in seq_entry['meta']:
                data_file_path = seq_entry['meta']['gvcf'].get('location')
                data_file_basename = seq_entry['meta']['gvcf'].get('basename')
                # Note this will handle only one secondary file.
                index_file_path = seq_entry['meta']['gvcf']['secondaryFiles'][0].get(
                    'location'
                )
                index_file_basename = seq_entry['meta']['gvcf']['secondaryFiles'][
                    0
                ].get('basename')

                sample_group_main = SampleGroup(
                    sample_id_external=external_sample_mapping.get(internal_sample_id),
                    sample_id_internal=internal_sample_id,
                    data_file=FileGroup(
                        path=data_file_path, basename=data_file_basename
                    ),
                    index_file=FileGroup(
                        path=index_file_path, basename=index_file_basename
                    ),
                    md5=FileGroup(
                        path=data_file_path + '.md5',
                        basename=data_file_basename + '.md5',
                    ),
                    batch_number=batch_number,
                )

                main_files.append(sample_group_main)

            if 'reads' in seq_entry['meta']:
                reads = seq_entry['meta']['reads']
                reads = reads if isinstance(reads, list) else [reads]
                for read in reads:
                    data_file_path = read.get('location')
                    data_file_basename = read.get('basename')
                    index_file_path = read['secondaryFiles'][0].get('location')
                    index_file_basename = read['secondaryFiles'][0].get('basename')
                    sample_group_archive = SampleGroup(
                        sample_id_external=external_sample_mapping.get(
                            internal_sample_id
                        ),
                        sample_id_internal=internal_sample_id,
                        data_file=FileGroup(
                            path=data_file_path, basename=data_file_basename
                        ),
                        index_file=FileGroup(
                            path=index_file_path, basename=index_file_basename
                        ),
                        md5=FileGroup(
                            path=data_file_path + '.md5',
                            basename=data_file_basename + '.md5',
                        ),
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


def create_analysis_in_sm_db(
    sample_group: SampleGroup, proj, output_path, analysis_type
):
    """ Creates a new analysis object"""
    aapi = AnalysisApi()

    internal_id = sample_group.sample_id_internal

    new_analysis = AnalysisModel(
        sample_ids=[internal_id],
        type=AnalysisType(analysis_type),
        status=AnalysisStatus('completed'),
        output=output_path,
    )

    aapi.create_new_analysis(proj, new_analysis)


def validate_md5(
    job: hb.batch.job, sample_group: SampleGroup, upload_path: str
) -> hb.batch.job:
    """Each gvcf and cram file is uploaded with a corresponding md5 checksum.
    This function appends commands to the job that moves these files, to validate each .md5 file."""

    base = sample_group.data_file.basename
    file_extension = '.'.join(base.split('.')[1:])
    file_path = f'{sample_group.sample_id_internal}.{file_extension}'

    # Generate paths to files that are being validated
    path_to_data = join(
        'gs://',
        upload_path,
        f'batch{sample_group.batch_number}',
        file_path,
    )

    path_to_md5 = join(
        'gs://',
        upload_path,
        f'batch{sample_group.batch_number}',
        f'{file_path}.md5',
    )

    # Calculate md5 checksum.
    job.command(
        f'gsutil cat {path_to_data} | md5sum | sed "s/-/{sample_group.data_file.basename}/" > /tmp/{sample_group.md5.basename}'
    )

    job.command(
        f'diff <(cat /tmp/{sample_group.md5.basename} | cut -d " " -f1 ) <(gsutil cat {path_to_md5} | cut -d " " -f1 )'
    )

    return job


def run_processor():
    """Set up and execute batch_move_files, validate_md5 and update_status"""

    # Setting up inputs for batch_move_files
    project = os.getenv('HAIL_BILLING_PROJECT')
    main_bucket = f'cpg-{project}-main'
    main_path = join(main_bucket, 'gvcf')
    archive_path = join(f'cpg-{project}-archive', 'cram')

    docker_image = os.environ.get('DRIVER_IMAGE')
    key = os.environ.get('GSA_KEY')

    if project is None:
        raise ValueError('HAIL_BILLING_PROJECT must be set')

    # Determine the analysis results (i.e. list of gvcfs and crams) to be moved
    main_files: List[SampleGroup] = []
    archive_files: List[SampleGroup] = []

    main_files, archive_files = determine_samples(project)

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
        output_path = os.path.join(
            'gs://',
            main_path,
            f'batch{sample_group.batch_number}',
            sample_group.sample_id_internal + '.g.vcf.gz',
        )

        status_job.call(
            create_analysis_in_sm_db,
            sample_group,
            project,
            output_path,
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
        output_path = os.path.join(
            'gs://',
            archive_path,
            f'batch{sample_group.batch_number}',
            sample_group.sample_id_internal + '.cram',
        )
        status_job.call(
            create_analysis_in_sm_db,
            sample_group,
            project,
            output_path,
            AnalysisType.CRAM,
        )
        status_job.depends_on(validate_job)
        archive_jobs.append(status_job)

    batch.run()


if __name__ == '__main__':
    run_processor()  # pylint: disable=E1120
