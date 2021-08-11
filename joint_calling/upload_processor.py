""" This function prepares gVCFs uploaded to GCS, based off sample status
    logged in a database, for further QC and downstream analysis.
    The upload processor will determine when samples should be added
    to existing MatrixTables where appropriate and which MatrixTables
    they should be combined with in this case. Following a successful
    run, all uploaded files will be moved to archival storage.

    Assumptions
    ===========
    - Specified files must exist in either the source or the destination
      bucket
    - Source and destination bucket must exist
    - User must be authenticated with appropriate permissions """

import os
from typing import List, Optional, NamedTuple
import hailtop.batch as hb


class FileGroup(NamedTuple):
    """Defines a file path and it's basename"""

    path: str
    basename: str


class SampleGroup(NamedTuple):
    """ Defines a group of files associated with each sample"""

    sample_id_external: str
    sample_id_internal: str
    data_file: FileGroup
    index_file: FileGroup
    md5: FileGroup
    batch_number: int


def batch_move_files(
    batch: hb.batch,
    sample_group: SampleGroup,
    destination_prefix: str,
    docker_image: Optional[str] = None,
    key: Optional[str] = None,
) -> List:
    """Creates a list of jobs to perform a batch move operation
    to move files from a source location to a destination.

    Parameters
    ==========
    batch: hb.Batch
        An object representing the DAG of jobs to run.
    files: SampleGroup
        A NamedTuple containing 3 files to be moved.
        For example ["TOB1543.g.vcf.gz","TOB1543.g.vcf.tbi","TOB1543.g.vcf.md5"]
        As well as the internal id, external id and batch number associated with each
        sample.
    source_prefix: str
        The path to the sub-directory where the files are initially located.
        For example "cpg-tob-wgs-upload" or "cpg-tob-wgs-upload/v1"
    destination_prefix: str
        The path to the sub-directory where the files should be moved.
        For example "cpg-tob-wgs-main" (not including the batch number)
    docker_image: str, optional
        The address and tag of a previously built docker image, within the
        artifact registry. Each batch job will run in this image.
        For example;
        australia-southeast1-docker.pkg.dev/project/images/driver:version'
    key: str, optional
        key-file for the service account used for authentication. In the case that this
        is not provided as an input, it is assumed that this key will exist at
        /gsa-key/key.json. This is the case when using the hail batch service backend.
        For example:
        "{
            "type": "service_account",
            "project_id": "",
            "private_key_id": "" ...
        }
        " """

    jobs = []

    external_id = sample_group.sample_id_external
    internal_id = sample_group.sample_id_internal
    batch_number = sample_group.batch_number

    files_to_move = (sample_group.data_file, sample_group.index_file, sample_group.md5)

    for file_name in files_to_move:
        base = file_name.basename
        previous_location = file_name.path

        # previous_location = os.path.join('gs://', source_prefix, file_name)
        # previous_location = file_name

        file_extension = base[len(external_id) :]
        new_file_name = internal_id + file_extension
        new_location = os.path.join(
            'gs://', destination_prefix, f'batch{batch_number}', new_file_name
        )

        j = batch.new_job(name=f'move {base} -> {new_file_name}')

        if docker_image is not None:
            j.image(docker_image)

        # Authenticate to service account.
        if key is not None:
            j.command(f"echo '{key}' > /tmp/key.json")
            j.command(
                f'gcloud -q auth activate-service-account --key-file=/tmp/key.json'
            )
        # Handles service backend, or a key in the same default location.
        else:
            j.command(
                'gcloud -q auth activate-service-account --key-file=/gsa-key/key.json'
            )

        # Checks file doesn't already exist at the destination, then performs move.
        j.command(
            f'gsutil -q stat "{new_location}" || '
            f'gsutil mv "{previous_location}" "{new_location}"'
        )
        jobs.append(j)

    return jobs
