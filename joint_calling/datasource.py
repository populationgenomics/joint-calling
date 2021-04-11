"""This module defines the methods required to authenticate,
access, filter data from a datasource (e.g. AirTable).
It will return a list of files to be moved by the upload
processor"""

import json
import abc
from typing import List
from airtable import Airtable
from google.cloud import secretmanager


class Datasource(metaclass=abc.ABCMeta):
    """A source of data"""

    @abc.abstractmethod
    def get_table(self):
        """Obtain a single table object from the
        data source"""

    @abc.abstractmethod
    def get_file_list(self, col, val) -> List[str]:
        """Returns a list of file names from the database based on a filter
        condition."""


@Datasource.register
class AirTableDatasource:
    """A series of functions implemented with AirTable as
    the datasource"""

    def __init__(self, project: str, secret_id: str):
        """Retrieves configuration file from secret manager
        in order to authenticate to airtable"""

        client = secretmanager.SecretManagerServiceClient()

        secret_name = f'projects/{project}/secrets/{secret_id}/versions/latest'
        response = client.access_secret_version(request={'name': secret_name})
        config_str = response.payload.data.decode('UTF-8')
        config = json.loads(config_str)
        self.project_config = config.get(project)

    def get_table(self) -> Airtable:
        """Authenticates to Airtable and returns the table specified
        in the config file."""

        # Get the Airtable credentials.
        base_key = self.project_config.get('baseKey')
        table_name = self.project_config.get('tableName')
        api_key = self.project_config.get('apiKey')

        # Pull table
        table = Airtable(base_key, table_name, api_key)

        return table

    def get_file_list(self, col: str, val: str, extension: str) -> List[str]:
        """Returns a list of file names from the database based on a filter
        condition."""
        table = self.get_table()
        sample_files = []
        for page in table.get_iter(formula=f"{col}='{val}'"):
            for record in page:
                sample = record['fields']['Sample_ID']
                sample_files.append(f'{sample}.{extension}')

        return sample_files
