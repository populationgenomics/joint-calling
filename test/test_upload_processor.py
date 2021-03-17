""" This file tests the functions defined in upload_processor """

import unittest
from unittest.mock import call, patch
from joint_calling.upload_processor import copy_files, delete_files


class TestUploadProcessor(unittest.TestCase):
    """A series of test cases for the upload processor"""

    def setUp(self):
        """ Initialises standard variables for testing """
        self.upload_bucket = 'test-dev-upload'
        self.main_bucket = 'test-dev-main'
        self.standard_list = ['Sample1.gVCF']
        self.invalid_bucket = 'invalid-bucket'

    @patch('joint_calling.upload_processor.storage')
    def test_copy(self, mock_storage):
        """Standard test case for the copy function, with standard input
        No errors should be produced nor exceptions raised. Passes in a standard
        list and assumes that the upload and main bucket are valid inputs"""

        copy_files(self.standard_list, self.upload_bucket, self.main_bucket)
        mock_storage.Client.assert_called_once()
        mock_client = mock_storage.Client()
        calls = [call.get_bucket(self.upload_bucket), call.get_bucket(self.main_bucket)]
        mock_client.assert_has_calls(calls, any_order=True)

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
