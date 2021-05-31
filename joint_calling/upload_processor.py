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
from typing import List, Optional
import hailtop.batch as hb


def batch_move_files(
    batch: hb.batch,
    files: List[str],
    source_prefix: str,
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
    files: List[str]
        A list of the file names to be moved.
        For example ["TOB1543.g.vcf.gz","TOB2314.g.vcf.gz","TOB3423.g.vcf.gz"]
    source_prefix: str
        The path to the sub-directory where the files are initially located.
        For example "cpg-tob-wgs-upload" or "cpg-tob-wgs-upload/v1"
    destination_prefix: str
        The path to the sub-directory where the files should be moved.
        For example "cpg-tob-wgs-main" or "cpg-tob-wgs-upload/batch1"
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
    for sample in files:
        previous_location = os.path.join('gs://', source_prefix, sample)
        new_location = os.path.join('gs://', destination_prefix, sample)

        j = batch.new_job(name=f'move {sample}')

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
