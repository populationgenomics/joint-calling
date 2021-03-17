""" This file tests the functions defined in upload_processor """

import unittest
from unittest.mock import call, patch
from google.api_core import exceptions
from joint_calling.upload_processor import copy_files, delete_files


class TestUploadProcessor(unittest.TestCase):
    """Test cases for the upload processor"""

    def setUp(self):
        """ Initialises standard variables for testing """
        self.upload_bucket = 'test-dev-upload'
        self.main_bucket = 'test-dev-main'
        self.standard_list = ['Sample1.gVCF']

    def test_copy(self):
        """A series of test cases for the copy function"""

        # Standard test case. Input list and buckets valid. Tests the get_buckets is
        # called twice, once for each input bucket.
        with patch('joint_calling.upload_processor.storage') as mock_storage:
            copy_files(self.standard_list, self.upload_bucket, self.main_bucket)
            mock_storage.Client.assert_called_once()
            mock_client = mock_storage.Client()
            calls = [
                call.get_bucket(self.upload_bucket),
                call.get_bucket(self.main_bucket),
            ]
            mock_client.assert_has_calls(calls, any_order=True)

        # Invalid bucket test case. Tests that an exception is raised when an invalid
        # bucket name is provided.
        with patch('joint_calling.upload_processor.storage') as mock_storage_invalid:
            mock_bucket = mock_storage_invalid.Client()
            mock_bucket.get_bucket.side_effect = exceptions.NotFound('Bucket Not Found')
            with self.assertRaises(Exception):
                copy_files(self.standard_list, self.upload_bucket, self.main_bucket)

    def test_delete(self):
        """A series of test cases for the delete function"""

        # Standard test case. Input list and buckets valid. Tests the appropriate
        # function calls
        with patch('joint_calling.upload_processor.storage') as mock_storage_delete:
            delete_files(self.standard_list, self.upload_bucket)
            mock_storage_delete.Client.assert_called_once()
            mock_client_delete = mock_storage_delete.Client()
            calls = [call.get_bucket(self.upload_bucket)]
            mock_client_delete.assert_has_calls(calls, any_order=True)

        # Invalid bucket test case. Tests that an exception is raised when an invalid
        # bucket name is provided.
        with patch(
            'joint_calling.upload_processor.storage'
        ) as mock_storage_invalid_delete:
            mock_bucket = mock_storage_invalid_delete.Client()
            mock_bucket.get_bucket.side_effect = exceptions.NotFound('Bucket Not Found')
            with self.assertRaises(Exception):
                delete_files(self.standard_list, self.upload_bucket)
