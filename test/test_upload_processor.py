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
        self.invalid_bucket = 'invalid-bucket'

    def test_copy(self):
        """A series of test cases for the copy function, using mock objects
        storage objects"""

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

    @patch('joint_calling.upload_processor.storage')
    def test_delete(self, mock_storage):
        """Standard test case for the copy function, with standard input
        No errors should be produced nor exceptions raised. Passes in a standard
        list and assumes that the upload and main bucket are valid inputs"""
        delete_files(self.standard_list, self.upload_bucket)
        mock_storage.Client.assert_called_once()
        mock_client = mock_storage.Client()
        calls = [call.get_bucket(self.upload_bucket)]
        mock_client.assert_has_calls(calls, any_order=True)
