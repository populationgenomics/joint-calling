"""A function to process uploaded samples"""

import json
import google.auth
from airtable import Airtable
from google.cloud import secretmanager


def get_table(project_config):
    """Get Table from Airtable"""

    # Get the Airtable credentials.
    base_key = project_config.get('baseKey')
    table_name = project_config.get('tableName')
    api_key = project_config.get('apiKey')

    # Pull table
    table = Airtable(base_key, table_name, api_key)

    return table


def get_config():
    """Retrieve configuration file from secret manager"""

    # secret_id set for testing, TODO: Replace
    secret_id = 'test-airtable-config'

    # Get project ID value. First output is a credentials object which isn't used.
    _, project = google.auth.default()

    client = secretmanager.SecretManagerServiceClient()

    secret_name = f'projects/{project}/secrets/{secret_id}/versions/latest'
    response = client.access_secret_version(request={'name': secret_name})
    config_str = response.payload.data.decode('UTF-8')
    config = json.loads(config_str)
    project_config = config.get(project)

    return project_config


if __name__ == '__main__':

    project_configuration = get_config()
    output_table = get_table(project_configuration)
