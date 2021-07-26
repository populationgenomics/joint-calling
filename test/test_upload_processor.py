#!/usr/bin/env python3

""" This file tests the functions defined in upload_processor """

import random
import string
import unittest
import subprocess
import os
from datetime import datetime
from typing import NamedTuple
import hailtop.batch as hb
from sample_metadata.api.sample_api import SampleApi
from joint_calling.upload_processor import batch_move_files, SampleGroup


class SuccessGroup(NamedTuple):
    """Similar to a SampleGroup, but excludes the index and md5 keys.
    Used to mimic a partial upload."""

    sample_id_external: str
    data_file: str


class FailedGroup(NamedTuple):
    """Similar to a SampleGroup, but excludes the index and md5 keys.
    Used to mimic a partial upload."""

    sample_id_external: str
    index_file: str
    md5: str


def validate_move(
    upload_prefix: str, main_prefix: str, original_file: str, new_file: str
) -> bool:
    """Checks that a given file exists at the location in
    the main bucket and no longer exists in the upload bucket.
    Returns True if this is the case and False otherwise.
    """

    main_path = os.path.join('gs://', main_prefix, new_file)
    upload_path = os.path.join('gs://', upload_prefix, original_file)
    exists_main = subprocess.run(['gsutil', '-q', 'stat', main_path], check=False)
    exists_upload = subprocess.run(['gsutil', '-q', 'stat', upload_path], check=False)

    # Exists at destination and not at source
    return exists_upload.returncode != 0 and exists_main.returncode == 0


def upload_files(sample_group: SampleGroup, upload_prefix: str):
    """Mimics the upload of a new sample. Takes a list of
    file names, creates these files, then moves them into
    the upload bucket used for further testing."""

    # Create and upload the files.

    for f in sample_group._fields:
        if f == 'sample_id_external':
            continue
        file_name = getattr(sample_group, f)
        subprocess.run(['touch', file_name], check=True)
        full_path = os.path.join('gs://', upload_prefix, file_name)
        subprocess.run(['gsutil', 'mv', file_name, full_path], check=True)


class TestUploadProcessor(unittest.TestCase):
    """Test cases for the upload processor"""

    def setUp(self):
        """Initialises standard variables for testing"""

        # Create random string of digits to label current test folder
        random_digits = ''.join(random.choice(string.digits) for i in range(5))
        timestamp = datetime.now().strftime('%d%m%Y_%H%M%S')
        test_folder = 'test' + random_digits + '_' + timestamp

        self.upload_prefix = os.path.join(
            'cpg-fewgenomes-test-tmp/test-upload', test_folder
        )
        self.main_prefix = os.path.join(
            'cpg-fewgenomes-test-tmp/test-main', test_folder
        )

        self.docker_image = os.environ.get('DOCKER_IMAGE')
        self.key = os.environ.get('GSA_KEY')
        self.project = 'viviandev'

    def test_batch_move_standard(self):
        """Testing standard case of moving a list of files with valid inputs"""

        # Assumption, that the following external ID already exists in the database.
        sapi = SampleApi()
        external_id = 'TOB02341'
        external_id_map = {'external_ids': [external_id]}
        internal_id_map = sapi.get_sample_id_map_by_external(
            self.project, external_id_map
        )
        internal_id = list(internal_id_map.values())[0]

        test_sample = SampleGroup(
            external_id,
            f'{external_id}.g.vcf.gz',
            f'{external_id}.g.vcf.gz.tbi',
            f'{external_id}.g.vcf.gz.md5',
        )

        upload_files(test_sample, self.upload_prefix)

        batch = hb.Batch(name='Test Batch Move Standard')
        batch_move_files(
            batch,
            test_sample,
            self.upload_prefix,
            self.main_prefix,
            self.project,
            self.docker_image,
            self.key,
        )

        batch.run()

        # Check that the files have been moved to main
        for f in test_sample._fields:
            if f == 'sample_id_external':
                continue
            file_name = getattr(test_sample, f)
            file_extension = file_name[len(test_sample.sample_id_external) :]
            self.assertTrue(
                validate_move(
                    self.upload_prefix,
                    self.main_prefix,
                    file_name,
                    str(internal_id) + file_extension,
                )
            )

    def test_batch_move_recovery(self):
        """Test cases that handles previous partially successful run.
        In this case, the file would not exist at the source
        but would exist at the destination"""

        # Rather than uploading files to the upload bucket, the sample
        # To the main bucket, to mimic this behaviour
        sapi = SampleApi()
        external_id = 'TOB02242'
        external_id_map = {'external_ids': [external_id]}
        internal_id_map = sapi.get_sample_id_map(self.project, external_id_map)
        internal_id = list(internal_id_map.values())[0]

        test_sample = SampleGroup(
            external_id,
            f'{external_id}.g.vcf.gz',
            f'{external_id}.g.vcf.gz.tbi',
            f'{external_id}.g.vcf.gz.md5',
        )

        recovery_sample = SampleGroup(
            str(internal_id),
            f'{internal_id}.g.vcf.gz',
            f'{internal_id}.g.vcf.gz.tbi',
            f'{internal_id}.g.vcf.gz.md5',
        )

        upload_files(recovery_sample, self.main_prefix)
        recovery_batch = hb.Batch(name='Test Batch Move Recovery')
        batch_move_files(
            recovery_batch,
            test_sample,
            self.upload_prefix,
            self.main_prefix,
            self.project,
            self.docker_image,
            self.key,
        )

        recovery_batch.run()

        # Check that the files have been moved to main
        for f in test_sample._fields:
            if f == 'sample_id_external':
                continue
            file_name = getattr(test_sample, f)
            file_extension = file_name[len(test_sample.sample_id_external) :]
            self.assertTrue(
                validate_move(
                    self.upload_prefix,
                    self.main_prefix,
                    file_name,
                    str(internal_id) + file_extension,
                )
            )

    def test_invalid_samples(self):
        """Test case that handles invalid sample ID's i.e. samples that don't exist
        in the upload bucket"""
        not_uploaded_id = 'TOB02337'
        not_uploaded_sample = SampleGroup(
            not_uploaded_id,
            f'{not_uploaded_id}.g.vcf.gz',
            f'{not_uploaded_id}.g.vcf.gz.tbi',
            f'{not_uploaded_id}.g.vcf.gz.md5',
        )
        invalid_batch = hb.Batch(name='Invalid Batch')
        batch_move_files(
            invalid_batch,
            not_uploaded_sample,
            self.upload_prefix,
            self.main_prefix,
            self.project,
            self.docker_image,
            self.key,
        )

        with self.assertRaises(subprocess.CalledProcessError):
            invalid_batch.run()

    def test_partial_recovery(self):
        """Another partial recovery test case. In this scenario,
        half the files were moved in a previous run and half were not."""

        sapi = SampleApi()
        external_id = 'TOB02341'
        external_id_map = {'external_ids': [external_id]}
        internal_id_map = sapi.get_sample_id_map(self.project, external_id_map)
        internal_id = list(internal_id_map.values())[0]

        full_sample = SampleGroup(
            external_id,
            f'{external_id}.g.vcf.gz',
            f'{external_id}.g.vcf.gz.tbi',
            f'{external_id}.g.vcf.gz.md5',
        )

        failed_files = FailedGroup(
            external_id,
            f'{external_id}.g.vcf.gz.tbi',
            f'{external_id}.g.vcf.gz.md5',
        )

        success_files = SuccessGroup(
            internal_id,
            f'{internal_id}.g.vcf.gz',
        )

        # Uploading this 'successful' files to the main bucket
        upload_files(success_files, self.main_prefix)

        # Uploading the 'failed' files to the upload bucket
        # i.e. replicating the case were they have not moved.
        upload_files(failed_files, self.upload_prefix)

        partial_batch = hb.Batch(name='Test Batch Move Partial Recovery')
        batch_move_files(
            partial_batch,
            full_sample,
            self.upload_prefix,
            self.main_prefix,
            self.project,
            self.docker_image,
            self.key,
        )

        partial_batch.run()

        # Check that the files have been moved to main
        for f in full_sample._fields:
            if f == 'sample_id_external':
                continue
            file_name = getattr(full_sample, f)
            file_extension = file_name[len(full_sample.sample_id_external) :]
            self.assertTrue(
                validate_move(
                    self.upload_prefix,
                    self.main_prefix,
                    file_name,
                    str(internal_id) + file_extension,
                )
            )

    def tearDown(self):
        """Deleting files created in test run"""
        full_path = os.path.join('gs://', self.main_prefix, '*')
        subprocess.run(['gsutil', 'rm', full_path], check=False)


if __name__ == '__main__':
    unittest.main()
