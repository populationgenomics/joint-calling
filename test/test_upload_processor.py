""" This file tests the functions defined in upload_processor """

import unittest
from unittest.mock import call, patch
from google.api_core import exceptions
from joint_calling.upload_processor import move_files


class TestUploadProcessor(unittest.TestCase):
    """Test cases for the upload processor"""

    def setUp(self):
        """ Initialises standard variables for testing """
        self.upload_bucket = 'test-dev-upload'
        self.main_bucket = 'test-dev-main'
        self.standard_list = ['Sample1.gVCF']
        self.standard_path = ''

    def test_move(self):
        """A series of test cases for the move function"""

        # Standard test case. Input list and buckets valid. Tests the get_buckets is
        # called twice, once for each input bucket.
        with patch('joint_calling.upload_processor.storage') as mock_storage:
            move_files(
                self.standard_list,
                self.upload_bucket,
                self.main_bucket,
                self.standard_path,
            )
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
                move_files(
                    self.standard_list,
                    self.upload_bucket,
                    self.main_bucket,
                    self.standard_path,
                )
