""" This file tests the functions defined in upload_processor """

import unittest
import subprocess
import os
from typing import List
import hailtop.batch as hb
from joint_calling.upload_processor import batch_move_files


def upload_files(files: List[str], upload_bucket: str):
    """A function to mimic file upload. Takes a list of
    file names, creates these files, then moved them into
    the upload bucket used for further testing."""

    for f in files:
        subprocess.run(['touch', f], check=True)
        subprocess.run(['gsutil', 'mv', f, f'gs://{upload_bucket}/{f}'], check=True)


class TestUploadProcessor(unittest.TestCase):
    """Test cases for the upload processor"""

    def setUp(self):
        """ Initialises standard variables for testing """
        self.upload_bucket = 'upload-processor-test-upload'
        self.main_bucket = 'upload-processor-test-main'
        self.batch_path = ''
        self.docker_image = os.environ.get('DOCKER_IMAGE')
        self.key = os.environ.get('KEY')

    def test_batch_move_standard(self):
        """ Testing standard case of moving a list of files with valid inputs"""
        sample_list = ['Sample5.gVCF', 'Sample6.gVCF', 'Sample7.gVCF']
        upload_files(sample_list, self.upload_bucket)
        batch = hb.Batch(name='Test Batch Move Standard')
        batch_move_files(
            batch,
            sample_list,
            self.upload_bucket,
            self.main_bucket,
            self.docker_image,
            self.key,
        )
        batch.run()
        # Check that the files have been moved to main
        for sample in sample_list:
            main_path = f'gs://{self.main_bucket}/{sample}'
            upload_path = f'gs://{self.upload_bucket}/{sample}'
            exists_main = subprocess.run(
                ['gsutil', '-q', 'stat', main_path], check=False
            )
            exists_upload = subprocess.run(
                ['gsutil', '-q', 'stat', upload_path], check=False
            )
            self.assertEqual(
                exists_main.returncode, 0, f'{sample} was not found in {main_path}'
            )
            self.assertEqual(
                exists_upload.returncode, 1, f'{sample} was found in {upload_path}'
            )

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
            main_path = f'gs://{self.main_bucket}/{sample}'
            upload_path = f'gs://{self.upload_bucket}/{sample}'
            exists_main = subprocess.run(
                ['gsutil', '-q', 'stat', main_path], check=False
            )
            exists_upload = subprocess.run(
                ['gsutil', '-q', 'stat', upload_path], check=False
            )
            self.assertEqual(
                exists_main.returncode, 0, f'{sample} was not found in {main_path}'
            )
            self.assertEqual(
                exists_upload.returncode, 1, f'{sample} was found in {upload_path}'
            )

    def test_invalid_samples(self):
        """Test case that handles invalid sample ID's i.e. samples that don't exist
        in the upload bucket"""
        sample_list = ['Sample8.gVCF', 'Sample9.gVCF']
        invalid_batch = hb.Batch(name='Invalid Batch')
        batch_move_files(
            invalid_batch,
            sample_list,
            self.upload_bucket,
            self.main_bucket,
            self.docker_image,
            self.key,
        )

        assertion_called = False

        try:
            invalid_batch.run()
        except subprocess.CalledProcessError:
            assertion_called = True

        self.assertTrue(assertion_called)
