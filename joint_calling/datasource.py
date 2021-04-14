"""This module defines the methods required to authenticate,
access, filter data from a datasource (e.g. AirTable).
It will return a list of files to be moved by the upload
processor"""

import json
import abc
from typing import List, Dict, Any
from enum import Enum
from airtable import Airtable
from google.cloud import secretmanager


class TableColumn(Enum):
    """Defines a series of table columns, and their corresponding ID's
    to be used for filtering"""

    STATUS = 'Status'
    SAMPLE_ID = 'Sample_ID'


class Datasource(metaclass=abc.ABCMeta):
    """A source of data"""

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

    def _get_table(self) -> Airtable:
        """Authenticates to Airtable and returns the table specified
        in the config file."""

        # Get the Airtable credentials.
        base_key = self.project_config.get('baseKey')
        table_name = self.project_config.get('tableName')
        api_key = self.project_config.get('apiKey')

        # Pull table
        table = Airtable(base_key, table_name, api_key)

        return table

    def get_file_list(
        self, filters: Dict[TableColumn, Any], extension: str
    ) -> List[str]:
        """Returns a list of file names from the database based on filter
        conditions.

        Parameters
        ==========
        filters: Dict[TableColumn, any]
        A dictionary where the key is a specific table column
        and the value represents the specific condition that needs to be matched.
        For example:
        {TableColumn.status: "Done"}
        This represents a filter on the status column, where status is set to
        Done.

        extension: str
        The file extension for the set of files to be returned.
        All files will be listed with this specific file extension.
        Assumes that this extension is valid & consistent for all samples.
        For example: '.gVCF', '.txt', etc.

        Returns
        =======


        """
        table = self._get_table()
        formula = get_formula(filters)
        sample_files = []
        for page in table.get_iter(formula=formula):
            for record in page:
                sample = record['fields']['Sample_ID']
                sample_files.append(f'{sample}.{extension}')

        return sample_files


def get_formula(filters: Dict[TableColumn, Any]) -> str:
    """Creates a forumula with the appropriate syntax for the datasource.
    Given a dictionary of filters, the function returns a formatted string
    in accordance with the guidelines of the datasource.
    """

    # In accordance with AirTable formula syntax
    formula = 'AND('

    # Append each filter and property to the formula
    for col_filter in filters:
        print(type(col_filter))
        formula = formula + col_filter.value + "='" + filters[col_filter] + "', "

    # Format & clean up
    formula = formula.strip(', ')
    formula = formula + ')'

    return formula
