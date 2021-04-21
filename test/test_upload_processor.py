""" This file tests the functions defined in upload_processor """

import unittest
import subprocess
import os
from typing import List
import hailtop.batch as hb
from joint_calling.upload_processor import batch_move_files


def validate_move(
    upload_bucket: str, main_bucket: str, sample: str, batch_path: str = ''
) -> bool:
    """Checks that a given sample file exists in the main bucket
    and no longer exists in the upload bucket

    Parameters
    ==========
    upload_bucket:
        str, the name of the upload bucket, for example "TOB-upload"
    main_bucket:
        str, the name of the main bucket, for example "fewgenomes-main"
    sample:
        the name of the file to be checked.
    batch_path:
        str, the path to the specific sub-directory where the file
        should be moved.


    Returns
    =======
    True when the file exists in the main bucket and not the upload bucket
    False otherwise.
    """

    upload_path = os.path.join('gs://', upload_bucket, sample)
    main_path = os.path.join('gs://', main_bucket, batch_path, sample)

    exists_main = subprocess.run(['gsutil', '-q', 'stat', main_path], check=False)
    exists_upload = subprocess.run(['gsutil', '-q', 'stat', upload_path], check=False)

    if exists_upload.returncode == 1 and exists_main.returncode == 0:
        return True

    return False


def cleanup(bucket: str, samples: List[str], path: str = ''):
    """Clean up function that deletes remaining files after
    test run"""

    for sample in samples:
        full_path = os.path.join('gs://', bucket, path, sample)
        subprocess.run(['gsutil', 'rm', full_path], check=False)


def upload_files(files: List[str], upload_bucket: str):
    """A function to mimic file upload. Takes a list of
    file names, creates these files, then moved them into
    the upload bucket used for further testing."""

    for f in files:
        subprocess.run(['touch', f], check=True)
        full_path = os.path.join('gs://', upload_bucket, f)
        subprocess.run(['gsutil', 'mv', f, full_path], check=True)


class TestUploadProcessor(unittest.TestCase):
    """Test cases for the upload processor"""

    def setUp(self):
        """ Initialises standard variables for testing """
        self.upload_bucket = 'upload-processor-test-upload'
        self.main_bucket = 'upload-processor-test-main'
        self.docker_image = os.environ.get('DOCKER_IMAGE')
        self.key = os.environ.get('KEY')

    def test_batch_move_standard(self):
        """ Testing standard case of moving a list of files with valid inputs"""
        sample_list = ['Sample5.gVCF', 'Sample6.gVCF', 'Sample7.gVCF']
        batch_path = 'batch0'
        upload_files(sample_list, self.upload_bucket)
        batch = hb.Batch(name='Test Batch Move Standard')
        batch_move_files(
            batch,
            sample_list,
            self.upload_bucket,
            self.main_bucket,
            self.docker_image,
            self.key,
            batch_path,
        )
        batch.run()
        # Check that the files have been moved to main
        for sample in sample_list:

            self.assertTrue(
                validate_move(self.upload_bucket, self.main_bucket, sample, batch_path)
            )

        cleanup(self.main_bucket, sample_list, batch_path)

    def test_batch_move_recovery(self):
        """Test cases that handles previous partially successful run.
        In this case, the file would not exist at the source
        but would exist at the destination"""

        # Rather than uploading files to the upload bucket, files are added
        # To the main bucket, to mimic this behaviour
        sample_list = ['Sample1.gVCF', 'Sample2.gVCF', 'Sample3.gVCF']
        upload_files(sample_list, self.main_bucket)
        recovery_batch = hb.Batch(name='Test Batch Move Recovery')
        batch_move_files(
            recovery_batch,
            sample_list,
            self.upload_bucket,
            self.main_bucket,
            self.docker_image,
            self.key,
        )
        recovery_batch.run()
        # Check that the files have been moved to main
        for sample in sample_list:
            self.assertTrue(validate_move(self.upload_bucket, self.main_bucket, sample))

        cleanup(self.main_bucket, sample_list)

    def test_invalid_samples(self):
        """Test case that handles invalid sample ID's i.e. samples that don't exist
        in the upload bucket"""
        sample_list = ['Sample8.gVCF', 'Sample9.gVCF']
        invalid_batch = hb.Batch(name='Invalid Batch')
        assertion_called = False
        try:
            batch_move_files(
                invalid_batch,
                sample_list,
                self.upload_bucket,
                self.main_bucket,
                self.docker_image,
                self.key,
            )
        except FileNotFoundError:
            assertion_called = True

        self.assertTrue(assertion_called)

    def test_partial_recovery(self):
        """Another partial recovery test case. In this scenario,
        half the files were moved in a previous run and half were not."""

        sample_list_failed = ['SampleA.gVCF', 'SampleB.gVCF']
        sample_list_successful = ['SampleC.gVCF', 'SampleD.gVCF']

        # Uploading this 'successful' files to the main bucket
        upload_files(sample_list_successful, self.main_bucket)

        # Uploading the 'failed' files to the upload bucket
        # i.e. replicating the case were they have not moved.
        upload_files(sample_list_failed, self.upload_bucket)

        sample_list = sample_list_failed + sample_list_successful

        partial_batch = hb.Batch(name='Test Batch Move Partial Recovery')
        batch_move_files(
            partial_batch,
            sample_list,
            self.upload_bucket,
            self.main_bucket,
            self.docker_image,
            self.key,
        )
        partial_batch.run()

        for sample in sample_list:
            self.assertTrue(validate_move(self.upload_bucket, self.main_bucket, sample))

        cleanup(self.main_bucket, sample_list)
