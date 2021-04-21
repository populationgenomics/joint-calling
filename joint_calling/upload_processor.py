""" This function prepares gVCF's uploaded to GCS, based off sample status
    logged in a database, for further QC and downstream analysis.
    The upload processor will determine when samples should be added
    to existing MatrixTables where appropriate and which MatrixTables
    they should be combined with in this case. Following a successful
    run all uploaded files will be moved to archival storage."""

import errno
import os
import subprocess
from typing import List
import hailtop.batch as hb


def batch_move_files(
    batch: hb.batch,
    sample_files: List[str],
    source_bucket: str,
    destination_bucket: str,
    docker_image: str,
    key: str = None,
    dest_path: str = '',
) -> List:
    """Moving files between buckets

    Parameters
    ==========
    batch: hb.Batch
        An object representing the DAG of jobs to run.
    sample_files: List[str]
        A list of the file names to be moved.
        For example ["TOB1543","TOB2314","TOB3423"]
    source_bucket: str
        The name of the bucket where the files are initially located.
        For example "cpg-tob-wgs-upload"
    destination_bucket: str
        The name of the bucket where files are to be moved.
        For example "cpg-tob-wgs-main"
    docker_image: str
        The address and tag of a previously built docker image, within the
        artifact registry.
        For example;
        australia-southeast1-docker.pkg.dev/project/images/driver:version'
    dest_path: str, optional
        The path to the specific sub-directory where the file should be moved.
        By default no sub-directory is specified. This may correspond to the version.
        For example, "/v1"
    key: str, optional
        key-file for the service account used for authentication. In the case that this
        is not provided as an input, it is assumed that this key will exist at
        gsa-key/key.json. This is the case when using the hail batch service backend.
        For example:
        "{
            "type": "service_account",
            "project_id": "",
            "private_key_id": "" ...
        }
        "

    Returns
    =======
    Returns a list of batch jobs. Each job consists of
    a gsutil move command for each valid file. When run,
    these jobs will perform the batch move file operation.

    Notes
    =====
    The function assumes that the file names have
    been validated at an earlier stage."""

    jobs = []

    for sample in sample_files:

        previous_location = os.path.join('gs://', source_bucket, sample)
        new_location = os.path.join('gs://', destination_bucket, dest_path, sample)

        # Checks if the files exist at the source and destination
        get_file_source = subprocess.run(
            ['gsutil', '-q', 'stat', previous_location], check=False
        )
        get_file_destination = subprocess.run(
            ['gsutil', '-q', 'stat', new_location], check=False
        )

        # Handles case where the file has already been moved in a previous run.
        # i.e. it exists at the destination and not the source.
        if get_file_source.returncode == 1 and get_file_destination.returncode == 0:
            continue

        # Invalid - the file does not exist at either location
        if get_file_source.returncode == 1 and get_file_destination.returncode == 1:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), sample)

        j = batch.new_job(name=f'moving_{sample}')
        j.image(docker_image)

        # Authenticate to service account.
        if key is not None:
            j.command(f"echo '{key}' > key.json")
            j.command(f'gcloud -q auth activate-service-account --key-file=key.json')
        # Handles service backend, or a key in the same default location.
        else:
            j.command(
                'gcloud -q auth activate-service-account --key-file=/gsa-key/key.json'
            )

        # -m performs a multi-threaded/multi-processing move
        j.command(f'gsutil -m mv {previous_location} {new_location}')
        jobs.append(j)

    return jobs
