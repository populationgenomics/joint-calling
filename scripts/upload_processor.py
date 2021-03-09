""" This function prepares gCVF's uploaded to GCP, based off sample status
    logged in AirTable, for joint calling. The upload processor will
    determine when samples should be added to existing MatrixTables where
    appropriate. Following a successful run all uploaded files will be
    moved to archival storage"""

import json
from airtable import Airtable
from google.cloud import secretmanager


def get_table(project_config: dict) -> Airtable:
    """Authenticates to Airtable and returns the table specified in the config file."""

    # Get the Airtable credentials.
    base_key = project_config.get('baseKey')
    table_name = project_config.get('tableName')
    api_key = project_config.get('apiKey')

    # Pull table
    table = Airtable(base_key, table_name, api_key)

    return table


def get_config(secret_id: str, project: str) -> dict:
    """Retrieves configuration file from secret manager when provided with a
    project and it's associated secret ID"""

    client = secretmanager.SecretManagerServiceClient()

    secret_name = f'projects/{project}/secrets/{secret_id}/versions/latest'
    response = client.access_secret_version(request={'name': secret_name})
    config_str = response.payload.data.decode('UTF-8')
    config = json.loads(config_str)
    project_config = config.get(project)

    return project_config
