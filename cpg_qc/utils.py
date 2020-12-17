import logging
import os
import subprocess
import sys
import time
from os.path import isdir, isfile, exists
from typing import Any, Callable

import click
from google.cloud import storage
import hashlib


def get_validation_callback(ext: str = None, must_exist: bool = False) -> Callable:
    def callback(ctx: click.Context, param: click.Option, value: Any):
        if value is None:
            return value
        if ext:
            assert isinstance(value, str), value
            value = value.rstrip('/')
            if not value.endswith(f'.{ext}'):
                raise click.BadParameter(
                    f'The argument {param.name} is expected to have '
                    f'an extension .{ext}, got: {value}')
        if must_exist:
            if not file_exists(value):
                raise click.BadParameter(
                    f'The value for {param.name} is expected to exist: {value}')
        return value
    return callback


def file_exists(path: str) -> bool:
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        path = path.replace('gs://', '').split('/', maxsplit=1)[1]
        if path.endswith('.mt') or path.endswith('.mt/'):
            path = os.path.join(path, '_SUCCESS')
        gs = storage.Client()
        return gs.get_bucket(bucket).get_blob(path)
    return os.path.exists(path)


def gs_cache_file(fpath: str, local_tmp_dir: str = None) -> str:
    """
    :param fpath: local or gs:// path. If the latter, the file will be downloaded
                  and cached if local_tmp_dir is provided, the local path will be
                  returned
    :param local_tmp_dir: a directory to cache files downloaded from cloud
    :return: file path
    """
    if fpath.startswith('gs://'):
        fname = os.path.basename(fpath) + '_' + hashlib.md5(fpath.encode()).hexdigest()[:6]
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
    """ Multiprocessing-safely and recursively creates a directory
    """
    if not dirpath:
        sys.stderr.write(f'Path is empty: {descriptive_name if descriptive_name else ""}\n')

    if isdir(dirpath):
        return dirpath

    if isfile(dirpath):
        sys.stderr.write(descriptive_name + ' ' + dirpath + ' is a file.\m')

    num_tries = 0
    max_tries = 10

    while not exists(dirpath):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dirpath)
        except OSError as e:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dirpath
