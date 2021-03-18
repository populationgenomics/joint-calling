""" This function prepares gCVF's uploaded to GCS, based off sample status
    logged in a database, for further QC and downstream analysis.
    The upload processor will determine when samples should be added
    to existing MatrixTables where appropriate and which MatrixTables
    they should be combined with in this case. Following a successful
    run all uploaded files will be moved to archival storage."""

from google.cloud import storage
from google.api_core.exceptions import NotFound


def copy_files(
    sample_files: list, source_bucket_name: str, target_bucket_name: str
) -> bool:
    """Given a list of file names, a source bucket and
    a target bucket, this function will copy the identified files
    to their new location. The function assumes that the file names have
    been validated at an earlier stage. It will return true when the
    files have been all successfully copied and false if at least 1 was not."""

    # Connecting to buckets
    try:
        storage_client = storage.Client()
        source_bucket = storage_client.get_bucket(source_bucket_name)
        target_bucket = storage_client.get_bucket(target_bucket_name)
    except NotFound as invalid_details:
        raise Exception('Bucket does not exist for this user.') from invalid_details

    current_files = storage_client.list_blobs(source_bucket)
    # A list of unique identifiers for each file in the bucket
    sample_crc32cs = []

    # Copying files from source bucket to target bucket
    for current_file in current_files:
        if current_file.name in sample_files:
            sample_crc32cs.append(current_file.crc32c)
            source_bucket.copy_blob(current_file, target_bucket)

    # Checks that the files were successfully copied over
    updated_files = storage_client.list_blobs(target_bucket)
    compare_list = []
    for updated_file in updated_files:
        compare_list.append(updated_file.crc32c)

    return set(sample_crc32cs).issubset(compare_list)


def delete_files(samples: list, target_bucket_name: str) -> bool:
    """Given a list of sample files and a target
    bucket, this function will delete the identified files
    from their specified location. This function assumes that
    the samples list is validated at an earlier stage. It will
    return true if all files are successfully deleted and false
    if at least one is not."""

    # Connecting to buckets

    try:
        storage_client = storage.Client()
        target_bucket = storage_client.get_bucket(target_bucket_name)
    except NotFound as invalid_details:
        raise Exception('Bucket does not exist for this user.') from invalid_details

    # A list of unique identifiers for each file in the bucket
    samples_crc32cs = []

    # Delete all files in the provided samples list
    current_files = storage_client.list_blobs(target_bucket)
    for current_file in current_files:
        if current_file.name in samples:
            samples_crc32cs.append(current_file.crc32c)
            target_bucket.delete_blob(current_file.name)

    # Checks that the deleted files no longer exist in the bucket
    deleted = True
    current_files = storage_client.list_blobs(target_bucket)
    for current_file in current_files:
        if current_file.crc32c in samples_crc32cs:
            deleted = False

    return deleted
