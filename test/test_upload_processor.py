import unittest
from joint_calling import upload_processor
from joint_calling.upload_processor import *
from unittest.mock import call, patch, MagicMock

class TestUploadProcessor(unittest.TestCase):
    def setUp(self):
        self.upload_bucket = "test-dev-upload"
        self.main_bucket = "test-dev-main"
        self.standard_list = ["Sample1.gVCF"]
        self.invalid_bucket = "invalid-bucket"

    @patch("joint_calling.upload_processor.storage")
    def test_copy(self, mock_storage):

        copy_files(self.standard_list, self.upload_bucket, self.main_bucket)
        mock_storage.Client.assert_called_once()
        mock_client = mock_storage.Client()
        calls = [call.get_bucket(self.upload_bucket),call.get_bucket(self.main_bucket)]
        mock_client.assert_has_calls(calls, any_order=True)

    @patch("joint_calling.upload_processor.storage")
    def test_delete(self, mock_storage):     

        delete_files(self.standard_list, self.upload_bucket)
        mock_storage.Client.assert_called_once()
        mock_client = mock_storage.Client()
        calls = [call.get_bucket(self.upload_bucket)]
        mock_client.assert_has_calls(calls, any_order=True)
