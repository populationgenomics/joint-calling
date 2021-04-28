#!/usr/bin/env python3

""" This file tests the functions defined in upload_processor """

import random
import string
import unittest
import subprocess
import os
from datetime import datetime
from typing import List
import hailtop.batch as hb
from joint_calling.upload_processor import batch_move_files


def validate_move(upload_prefix: str, main_prefix: str, sample: str) -> bool:
    """Checks that a given file exists at the location in
    the main bucket and no longer exists in the upload bucket.
    Returns True if this is the case and False otherwise.
    """

    main_path = os.path.join('gs://', main_prefix, sample)
    upload_path = os.path.join('gs://', upload_prefix, sample)

    exists_main = subprocess.run(['gsutil', '-q', 'stat', main_path], check=False)
    exists_upload = subprocess.run(['gsutil', '-q', 'stat', upload_path], check=False)

    # Exists at destination and not at source
    return exists_upload.returncode != 0 and exists_main.returncode == 0


def upload_files(files: List[str], upload_prefix: str):
    """A function to mimic file upload. Takes a list of
    file names, creates these files, then moves them into
    the upload bucket used for further testing."""

    for f in files:
        subprocess.run(['touch', f], check=True)
        full_path = os.path.join('gs://', upload_prefix, f)
        subprocess.run(['gsutil', 'mv', f, full_path], check=True)


class TestUploadProcessor(unittest.TestCase):
    """Test cases for the upload processor"""

    def setUp(self):
        """Initialises standard variables for testing"""

        # Create random string of digits to label current test folder
        random_digits = ''.join(random.choice(string.digits) for i in range(5))
        timestamp = datetime.now().strftime('%d%m%Y_%H%M%S')
        test_folder = 'test' + random_digits + '_' + timestamp

        self.upload_prefix = os.path.join(
            'cpg-fewgenomes-temporary/test-upload', test_folder
        )
        self.main_prefix = os.path.join(
            'cpg-fewgenomes-temporary/test-main', test_folder
        )
        self.docker_image = os.environ.get('DOCKER_IMAGE')
        self.key = os.environ.get('GSA_KEY')

    def test_batch_move_standard(self):
        """Testing standard case of moving a list of files with valid inputs"""
        sample_list = ['Sample5.gVCF', 'Sample6.gVCF', 'Sample7.gVCF']
        upload_files(sample_list, self.upload_prefix)
        batch = hb.Batch(name='Test Batch Move Standard')
        batch_move_files(
            batch,
            sample_list,
            self.upload_prefix,
            self.main_prefix,
            self.docker_image,
            self.key,
        )

        batch.run()

        # Check that the files have been moved to main
        for sample in sample_list:
            self.assertTrue(validate_move(self.upload_prefix, self.main_prefix, sample))

    def test_batch_move_recovery(self):
        """Test cases that handles previous partially successful run.
        In this case, the file would not exist at the source
        but would exist at the destination"""

        # Rather than uploading files to the upload bucket, files are added
        # To the main bucket, to mimic this behaviour
        sample_list = ['Sample1.gVCF', 'Sample2.gVCF', 'Sample3.gVCF']
        upload_files(sample_list, self.main_prefix)
        recovery_batch = hb.Batch(name='Test Batch Move Recovery')
        batch_move_files(
            recovery_batch,
            sample_list,
            self.upload_prefix,
            self.main_prefix,
            self.docker_image,
            self.key,
        )

        recovery_batch.run()

        # Check that the files have been moved to main
        for sample in sample_list:
            self.assertTrue(validate_move(self.upload_prefix, self.main_prefix, sample))

    def test_invalid_samples(self):
        """Test case that handles invalid sample ID's i.e. samples that don't exist
        in the upload bucket"""
        sample_list = ['Sample8.gVCF', 'Sample9.gVCF']
        invalid_batch = hb.Batch(name='Invalid Batch')
        batch_move_files(
            invalid_batch,
            sample_list,
            self.upload_prefix,
            self.main_prefix,
            self.docker_image,
            self.key,
        )

        with self.assertRaises(subprocess.CalledProcessError):
            invalid_batch.run()

    def test_partial_recovery(self):
        """Another partial recovery test case. In this scenario,
        half the files were moved in a previous run and half were not."""

        sample_list_failed = ['SampleA.gVCF', 'SampleB.gVCF']
        sample_list_successful = ['SampleC.gVCF', 'SampleD.gVCF']

        # Uploading this 'successful' files to the main bucket
        upload_files(sample_list_successful, self.main_prefix)

        # Uploading the 'failed' files to the upload bucket
        # i.e. replicating the case were they have not moved.
        upload_files(sample_list_failed, self.upload_prefix)

        sample_list = sample_list_failed + sample_list_successful

        partial_batch = hb.Batch(name='Test Batch Move Partial Recovery')
        batch_move_files(
            partial_batch,
            sample_list,
            self.upload_prefix,
            self.main_prefix,
            self.docker_image,
            self.key,
        )

        partial_batch.run()

        for sample in sample_list:
            self.assertTrue(validate_move(self.upload_prefix, self.main_prefix, sample))

    def tearDown(self):
        # Deleting files created in test run
        full_path = os.path.join('gs://', self.main_prefix, '*')
        subprocess.run(['gsutil', 'rm', full_path], check=False)


if __name__ == '__main__':
    unittest.main()
