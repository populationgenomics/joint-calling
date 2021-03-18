""" This function prepares gVCF's uploaded to GCS, based off sample status
    logged in a database, for further QC and downstream analysis.
    The upload processor will determine when samples should be added
    to existing MatrixTables where appropriate and which MatrixTables
    they should be combined with in this case. Following a successful
    run all uploaded files will be moved to archival storage."""

from typing import Tuple, List
from google.cloud import storage
from google.api_core.exceptions import NotFound


def copy_files(
    sample_files: List[str], source_bucket_name: str, destination_bucket_name: str
) -> Tuple[List[str], List[str]]:

    """Given a list of file names, a source bucket and
    a destination bucket, this function will copy the identified files
    to their new location. The function assumes that the file names have
    been validated at an earlier stage. It will return two lists. One containing
    the successfully copied files, and one containing the files that were not copied."""

    # Connecting to buckets
    try:
        storage_client = storage.Client()
        source_bucket = storage_client.get_bucket(source_bucket_name)
        destination_bucket = storage_client.get_bucket(destination_bucket_name)
    except NotFound as invalid_details:
        raise Exception('Bucket does not exist for this user.') from invalid_details

    skipped = []
    copied = []

    # Copying files from source to destination
    for sample in sample_files:
        sample_blob = source_bucket.get_blob(sample)
        if sample_blob:
            # Check if the sample also exists in the destination folder.
            check_blob = destination_bucket.get_blob(sample)
            if check_blob:
                if check_blob.crc32c == sample_blob.crc32c:
                    # skip the copy as the same file already exists in the
                    # destination folder
                    skipped.append(sample)
                else:
                    # A file with the same name exists, however it contains different
                    # content, so the new version is copied across
                    source_bucket.copy_blob(sample_blob, destination_bucket)
                    copied.append(sample)
            else:
                # copy the file over.
                source_bucket.copy_blob(sample_blob, destination_bucket)
                copied.append(sample)
        else:
            # Sample file was not found
            skipped.append(sample)

    return copied, skipped


def delete_files(samples: List[str], bucket_name: str) -> Tuple[List[str], List[str]]:
    """Given a list of sample files and a
    bucket, this function will delete the identified files
    from their specified location. This function assumes that
    the samples list is validated at an earlier stage. It will return
    two lists. One containing the successfully copied files, and one
    containing the unsuccessfully copied files."""

    # Connecting to buckets
    try:
        storage_client = storage.Client()
        target_bucket = storage_client.get_bucket(bucket_name)
    except NotFound as invalid_details:
        raise Exception('Bucket does not exist for this user.') from invalid_details

    skipped = []
    deleted = []

    for sample in samples:
        try:
            # Delete the specified sample
            target_bucket.delete_blob(sample)
            deleted.append(sample)
        except NotFound:
            # The sample file was not found
            skipped.append(sample)

    return deleted, skipped
