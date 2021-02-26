"""Utility functions for the cpg-qc module"""

import os
import sys
import time
import hashlib
from os.path import isdir, isfile, exists
from typing import Any, Callable

import hail as hl
import click
from google.cloud import storage


DEFAULT_REF = 'GRCh38'


def init_hail(name: str, local_tmp_dir: str):
    """
    Initialize Hail and set up the directory for logs
    :param name: name to prefix the log file
    :param local_tmp_dir: local directory to write Hail logs
    :return:
    """
    timestamp = time.strftime('%Y%m%d-%H%M')
    hl_log = os.path.join(
        safe_mkdir(os.path.join(local_tmp_dir, 'log')), f'{name}-{timestamp}.log'
    )
    hl.init(default_reference=DEFAULT_REF, log=hl_log)


def init_hail_service(
    name: str, work_bucket: str, hail_billing_project: str, local_tmp_dir: str
):
    """
    Initialize Hail and set up the directory for logs
    :param name: name to prefix the log file
    :param local_tmp_dir: local directory to write Hail logs
    :return:
    """
    timestamp = time.strftime('%Y%m%d-%H%M')
    hl_log = os.path.join(
        safe_mkdir(os.path.join(local_tmp_dir, 'log')), f'{name}-{timestamp}.log'
    )
    hl.init(default_reference=DEFAULT_REF, log=hl_log)
    hl.context.init_service(
        default_reference=DEFAULT_REF,
        bucket=work_bucket,
        billing_project=hail_billing_project,
    )


def get_validation_callback(
    ext: str = None,
    must_exist: bool = False,
    accompanying_metadata_suffix: str = None,
) -> Callable:
    """
    Get callback for Click parameters validation
    :param ext: check that the path has the expected extention
    :param must_exist: check that the input file/object/directory exists
    :param accompanying_metadata_suffix: checks that a file at the same location but
    with a different suffix also exists (e.g. genomes.mt and genomes.metadata.ht)
    :return: a callback suitable for Click parameter initialization
    """

    def callback(_: click.Context, param: click.Option, value: Any):
        if value is None:
            return value
        if ext:
            assert isinstance(value, str), value
            value = value.rstrip('/')
            if not value.endswith(f'.{ext}'):
                raise click.BadParameter(
                    f'The argument {param.name} is expected to have '
                    f'an extension .{ext}, got: {value}'
                )
        if must_exist:
            if not file_exists(value):
                raise click.BadParameter(f"{value} doesn't exist or incomplete")
            if accompanying_metadata_suffix:
                accompanying_metadata_fpath = (
                    os.path.splitext(value)[0] + accompanying_metadata_suffix
                )
                if not file_exists(accompanying_metadata_fpath):
                    raise click.BadParameter(
                        f"An accompanying file {accompanying_metadata_fpath} doesn't "
                        f'exist'
                    )
        return value

    return callback


def file_exists(path: str) -> bool:
    """
    Check if the object exists, where the object can be:
        * local file
        * local directory
        * Google Storage object
        * Google Storage URL representing a *.mt or *.ht Hail data,
          in which case it will check for the existince of a
          *.mt/_SUCCESS or *.ht/_SUCCESS file.
    :param path: path to the file/directory/object/mt/ht
    :return: True if the object exists
    """
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        path = path.replace('gs://', '').split('/', maxsplit=1)[1]
        path = path.rstrip('/')  # ".mt/" -> ".mt"
        if any(path.endswith(f'.{suf}') for suf in ['mt', 'ht']):
            path = os.path.join(path, '_SUCCESS')
        gs = storage.Client()
        return gs.get_bucket(bucket).get_blob(path)
    return os.path.exists(path)


def gs_cache_file(fpath: str, local_tmp_dir: str = None) -> str:
    """
    :param fpath: local or a `gs://` path. If the latter, the file
        will be downloaded and cached if local_tmp_dir is provided,
        the local path will be returned
    :param local_tmp_dir: a local directory to cache files downloaded
        from Google Storage
    :return: file path
    """
    if fpath.startswith('gs://'):
        fname = (
            os.path.basename(fpath) + '_' + hashlib.md5(fpath.encode()).hexdigest()[:6]
        )
        local_fpath = os.path.join(local_tmp_dir, fname)
        if not exists(local_fpath):
            bucket = fpath.replace('gs://', '').split('/')[0]
            path = fpath.replace('gs://', '').split('/', maxsplit=1)[1]
            gs = storage.Client()
            blob = gs.get_bucket(bucket).get_blob(path)
            if blob:
                blob.download_to_filename(local_fpath)
    else:
        local_fpath = fpath
    return local_fpath


def safe_mkdir(dirpath: str, descriptive_name: str = '') -> str:
    """
    Multiprocessing-safely and recursively creates a directory
    """
    if not dirpath:
        sys.stderr.write(
            f'Path is empty: {descriptive_name if descriptive_name else ""}\n'
        )

    if isdir(dirpath):
        return dirpath

    if isfile(dirpath):
        sys.stderr.write(descriptive_name + ' ' + dirpath + ' is a file.\n')

    num_tries = 0
    max_tries = 10

    while not exists(dirpath):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dirpath)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dirpath
