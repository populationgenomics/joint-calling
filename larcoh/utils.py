"""
Utility functions and constants.
"""
import collections
import os
import re
import unicodedata
import logging
import string
import sys
import time
import traceback
from functools import lru_cache
from random import choices
from typing import cast, Callable

import hail as hl
import click

from cpg_utils import Path, to_path
from cpg_utils.config import get_config

logger = logging.getLogger(__file__)


# Packages to install on a dataproc cluster, to use with the dataproc wrapper.
DATAPROC_PACKAGES = [
    'cpg_utils',
    'larcoh',
    'cpg_gnomad',  # github.com/populationgenomics/gnomad_methods
    'seqr_loader==1.2.5',  # hail-elasticsearch-pipelines
    'elasticsearch==8.1.1',
    'cpg_utils',
    'click',
    'google',
    'fsspec',
    'sklearn',
    'gcloud',
    'selenium',
]


@lru_cache
def exists(
    path: Path | str, verbose: bool = True, description: str | None = None
) -> bool:
    """
    Caching version of the existence check.
    The python code runtime happens entirely during the pipeline submission,
    without waiting for it to finish, so there is no expectation that object
    existence status would change during the runtime. This, this function uses
    `@lru_cache` to make sure that object existence is checked only once.
    """
    return exists_not_cached(path, verbose)


def exists_not_cached(
    path: Path | str, verbose: bool = True, description: str | None = None
) -> bool:
    """
    Check if the object exists, where the object can be:
        * local file
        * local directory
        * cloud object
        * cloud URL representing a *.mt or *.ht Hail data,
          in which case it will check for the existence of a
          *.mt/_SUCCESS or *.ht/_SUCCESS file.
    @param path: path to the file/directory/object/mt/ht
    @param verbose: print on each check
    @param description: optional string to print in verbose mode
    @return: True if the object exists
    """
    extra_msg = '' if not description else f' ({description})'

    path = cast(Path, to_path(path))

    # rstrip to ".mt/" -> ".mt"
    if any(str(path).rstrip('/').endswith(f'.{suf}') for suf in ['mt', 'ht']):
        path = path / '_SUCCESS'

    if verbose:
        # noinspection PyBroadException
        try:
            res = path.exists()
        except BaseException:
            traceback.print_exc()
            logger.error(f'Failed checking {path}' + extra_msg)
            sys.exit(1)
        logger.debug(
            f'Checked {path} [' + ('exists' if res else 'missing') + ']' + extra_msg
        )
        return res
    return path.exists()


def can_reuse(path: list[Path] | Path | str | None) -> bool:
    """
    Checks if `fpath` is good to reuse in the analysis: it exists
    and `overwrite` is False.

    If `fpath` is a collection, it requires all files in it to exist.
    """
    if get_config()['workflow'].get('overwrite', False):
        return False

    if not path:
        return False

    if isinstance(path, list):
        return all(can_reuse(fp) for fp in path)

    if not exists(path):
        return False

    logger.debug(f'Reusing existing {path}. Use --overwrite to overwrite')
    return True


def timestamp(rand_suffix_len: int = 5) -> str:
    """
    Generate a timestamp string. If `rand_suffix_len` is set, adds a short random
    string of this length for uniqueness.
    """
    result = time.strftime('%Y_%m%d_%H%M')
    if rand_suffix_len:
        rand_bit = ''.join(
            choices(string.ascii_uppercase + string.digits, k=rand_suffix_len)
        )
        result += f'_{rand_bit}'
    return result


def slugify(line: str):
    """
    Slugify a string.

    Example:
    >>> slugify(u"Héllø Wörld")
    u"hello-world"
    """

    return re.sub(
        r'[-\s]+',
        '-',
        re.sub(
            r'[^\w\s-]',
            '',
            unicodedata.normalize('NFKD', line).encode('ascii', 'ignore').decode(),
        )
        .strip()
        .lower(),
    )


def check_duplicates(iterable):
    """
    Throws error if input list contains repeated items.
    """
    duplicates = [
        item for item, count in collections.Counter(iterable).items() if count > 1
    ]
    if duplicates:
        raise ValueError(f'Found {len(duplicates)} duplicates: {duplicates}')
    return duplicates


def get_validation_callback(
    ext: str = None,
    must_exist: bool = False,
    accompanying_metadata_suffix: str = None,
) -> Callable:
    """
    Get callback for Click parameters validation
    :param ext: check that the path has the expected extension
    :param must_exist: check that the input file/object/directory exists
    :param accompanying_metadata_suffix: checks that a file at the same location but
    with a different suffix also exists (e.g. genomes.mt and genomes.metadata.ht)
    :return: a callback suitable for Click parameter initialization
    """

    def callback(_, param, value):
        if value is None:
            return None
        if ext:
            assert isinstance(value, str), value
            value = value.rstrip('/')
            if not value.endswith(f'.{ext}'):
                raise click.BadParameter(
                    f'The argument {param.name} is expected to have '
                    f'an extension .{ext}, got: {value}'
                )
        if must_exist:
            if not exists(value):
                raise click.BadParameter(f"{value} doesn't exist or incomplete")
            if accompanying_metadata_suffix:
                accompanying_metadata_fpath = (
                    os.path.splitext(value)[0] + accompanying_metadata_suffix
                )
                if not exists(accompanying_metadata_fpath):
                    raise click.BadParameter(
                        f"An accompanying file {accompanying_metadata_fpath} doesn't "
                        f'exist'
                    )
        return value

    return callback


def get_mt(
    mt_path: str,
    split: bool = False,
    hard_filtered_samples_to_remove_ht: hl.Table = None,
    meta_ht: hl.Table = None,
    add_meta: bool = False,
    release_only: bool = False,
    passing_sites_only: bool = False,
    unrelated_only: bool = False,
) -> hl.MatrixTable:
    """
    Wrapper function to get data with desired filtering and metadata annotations
    :param mt_path: path to the MatrixTable
    :param split:
        Split multiallelics and convert local-allele LGT/LA fields to GT.
        Note: this will perform a split on the MT rather than grab an already split MT
    :param hard_filtered_samples_to_remove_ht:
        table with samples to remove
        (only relevant after sample QC that produces a table with samples failed
        filtering)
    :param meta_ht: table with meta-information generated by sample QC
    :param add_meta: whether to add metadata to MT in 'meta' column
    :param release_only: whether to filter the MT to only samples available for
        release (can only be used with)
    :param passing_sites_only: whether to filter the MT to only variants with
        nothing in the filter field (e.g. passing soft filters)
    :param unrelated_only: remove related samples (keep one sample from a family)
    :return: MatrixTable with chosen annotations and filters
    """
    mt = hl.read_matrix_table(mt_path)

    if passing_sites_only:
        try:
            mt = mt.filter_rows(mt.filters.length() == 0)
        except AttributeError:
            pass

    if hard_filtered_samples_to_remove_ht is not None:
        mt = mt.filter_cols(
            hl.is_missing(hard_filtered_samples_to_remove_ht[mt.col_key])
        )

    if add_meta:
        assert meta_ht is not None
        mt = mt.annotate_cols(meta=meta_ht[mt.col_key])

        if release_only:
            mt = mt.filter_cols(mt.meta.release)

        if unrelated_only:
            mt = mt.filter_cols(~mt.meta.related)

    else:
        if release_only:
            assert meta_ht is not None
            mt = mt.filter_cols(meta_ht[mt.col_key].release)

        if unrelated_only:
            assert meta_ht is not None
            mt = mt.filter_cols(~meta_ht[mt.col_key].related)

    if split:
        mt = mt.annotate_rows(
            n_unsplit_alleles=hl.len(mt.alleles),
            mixed_site=(hl.len(mt.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(mt.alleles[0], a), mt.alleles[1:])
            & hl.any(lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:]),
        )
        # Will use GT instead of LGT
        mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    return mt
