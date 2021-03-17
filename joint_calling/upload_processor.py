""" This function prepares gCVF's uploaded to GCP, based off sample status
    logged in a database, for joint calling. The upload processor will
    determine when samples should be added to existing MatrixTables where
    appropriate and which MatrixTables they should be combined with.
    Following a successful run all uploaded files will be moved to archival 
    storage"""

import json
from typing import Tuple
from google.cloud import storage

def copy_files(sample_files: list, source_bucket_name: str, target_bucket_name: str) -> bool:
    """ Given a list of file names, an source bucket and 
    a target bucket, this function will copy the identified files
    to their new location. Assumes that sample_files has been validated 
    at an earlier stage. Assumes source bucket has been validated."""

    #Connecting to buckets.
    try: 
        storage_client = storage.Client()
        source_bucket = storage_client.get_bucket(source_bucket_name)
        target_bucket = storage_client.get_bucket(target_bucket_name)
    except:
        raise Exception("Invalid GCS details. Double check your bucket names & project config")
        return False

    
    #Copying files from source bucket to target bucket
    current_files = storage_client.list_blobs(source_bucket)
    sample_crc32cs = []

    for current_file in current_files:
        if(current_file.name in sample_files):
            sample_crc32cs.append(current_file.crc32c)
            source_bucket.copy_blob(current_file, target_bucket)

    # Checks that the files were successfully copied over. 
    updated_files = storage_client.list_blobs(target_bucket)
    compare_list = []
    for updated_file in updated_files:
        compare_list.append(updated_file.crc32c)
    
    return set(sample_crc32cs).issubset(compare_list)


def delete_files(samples:list, target_bucket_name: str) -> bool:
    """ Given a list of sample files and a target 
    bucket, this function will delete the identified files
    from their specified location """
    #Connecting to buckets.

    deleted = True

    try: 
        storage_client = storage.Client()
        target_bucket = storage_client.get_bucket(target_bucket_name)
    except:
        raise Exception("Invalid GCS details. Double check your bucket names & project config")
        return False
    
    #Delete all files in the provided samples list
    current_files = storage_client.list_blobs(target_bucket)
    for current_file in current_files:
        if(current_file.name in samples):
            target_bucket.delete_blob(current_file.name)
    
    #Check that they were deleted
    current_files = storage_client.list_blobs(target_bucket)
    for current_file in current_files:
        if(current_file.name in samples):
            deleted = False

    return deleted
