""" This function prepares gVCF's uploaded to GCS, based off sample status
    logged in a database, for further QC and downstream analysis.
    The upload processor will determine when samples should be added
    to existing MatrixTables where appropriate and which MatrixTables
    they should be combined with in this case. Following a successful
    run all uploaded files will be moved to archival storage."""

from typing import List
from google.cloud import storage
from google.api_core.exceptions import NotFound


def move_files(
    sample_files: List[str],
    source_bucket_name: str,
    destination_bucket_name: str,
    path: str = '',
):
    """Moving files between buckets

    Parameters
    ==========
    sample_files: List[str]
        A list of the file names to be moved.
        For example ["TOB1543","TOB2314","TOB3423"]
    source_bucket_name: str
        The name of the bucket where the files are initially located.
        For example "cpg-tob-wgs-upload"
    destination_bucket_name: str
        The name of the bucket where files are to be moved.
        For example "cpg-tob-wgs-main"
    path: str, optional
        The path to the specific sub-directory where the file should be moved.
        By default no sub-directory is specified.
        For example, "/v1"

    Notes
    =====
    The function assumes that the file names have
    been validated at an earlier stage.

    The function will raise an error if invalid bucket details,
    or invalid sample names are produced."""

    # Connecting to buckets
    try:
        storage_client = storage.Client()
        source_bucket = storage_client.get_bucket(source_bucket_name)
        destination_bucket = storage_client.get_bucket(destination_bucket_name)
    except NotFound as invalid_details:
        print(invalid_details)
        raise NotFound('Invalid bucket details') from invalid_details

    # Iterating through files to be moved
    for sample in sample_files:
        sample_blob = source_bucket.get_blob(sample)
        destination_blob = destination_bucket.get_blob(path + sample)
        # Sample exists in both the source and the destination buckets.
        if destination_blob and sample_blob:
            # Check if the files have identical content
            if destination_blob.crc32c == sample_blob.crc32c:
                # skip the copy as the same file already exists in the
                # destination folder
                source_bucket.delete_blob(sample_blob.name)
            else:
                # A file with the same name exists, however it contains different
                # content, so the new version is moved across
                destination_bucket.rename_blob(sample_blob, path + sample_blob.name)
        # Sample file has not been added to the destination previously
        elif sample_blob and not destination_blob:
            # Standard case - move the file
            destination_bucket.rename_blob(sample_blob, path + sample_blob.name)
        elif destination_blob and not sample_blob:
            # The sample was moved to the destination in a previous execution
            # of the pipeline.
            continue
        else:
            raise NotFound(f'{sample} was not found')
