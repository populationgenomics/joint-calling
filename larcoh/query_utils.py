"""
Utility functions for Query script submitted with Dataproc. Do not depend on
cpg-pipes nad python 3.10.
"""

import collections
import logging
import os
import sys
import traceback
from functools import lru_cache
from typing import Callable, List, Union, Optional, cast

import click
import hail as hl
from cpg_utils import Path, to_path
from cpg_utils.config import get_config

logger = logging.getLogger(__file__)


@lru_cache
def exists(path: Union[Path, str], verbose: bool = True) -> bool:
    """
    Caching version of the existence check.
    The python code runtime happens entirely during the pipeline submission,
    without waiting for it to finish, so there is no expectation that object
    existence status would change during the runtime. This, this function uses
    `@lru_cache` to make sure that object existence is checked only once.
    """
    return exists_not_cached(path, verbose)


def exists_not_cached(path: Optional[Path], verbose: bool = True) -> bool:
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
    @return: True if the object exists
    """
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
            logger.error(f'Failed checking {path}')
            sys.exit(1)
        logger.debug(f'Checked {path} [' + ('exists' if res else 'missing') + ']')
        return res
    return path.exists()


def can_reuse(path: Union[List[Path], Path, str, None]) -> bool:
    """
    Checks if `path` is good to reuse in the analysis: it exists,
    and `overwrite` is False.

    If `fpath` is a collection, it requires all files in it to exist.
    """
    overwrite = get_config()['workflow'].get('overwrite', False)
    if overwrite is True:
        return False

    if not path:
        return False

    if isinstance(path, list):
        return all(can_reuse(fp) for fp in path)

    if not exists(path):
        return False

    logger.debug(f'Reusing existing {path}. Use --overwrite to overwrite')
    return True


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


def get_vds(
    vds_path: str,
    split: bool = False,
    hard_filtered_samples_to_remove_ht: Optional[hl.Table] = None,
    meta_ht: Optional[hl.Table] = None,
    release_only: bool = False,
    n_partitions: int = None,
) -> hl.vds.VariantDataset:
    """
    Wrapper function to get VDS with desired filtering and metadata annotations.
    @param vds_path: path to VDS
    @param split:
        Split multiallelics and convert local-allele LGT/LA fields to GT.
        Note: this will perform a split on the MT rather than grab an already split MT
    @param hard_filtered_samples_to_remove_ht:
        table with samples to remove
        (only relevant after sample QC that produces a table with samples failed
        filtering)
    @param meta_ht: table with meta-information generated by sample QC
    @param release_only: whether to filter the MT to only samples available for
        release (can only be used with)
    @param n_partitions: number of partitions to read the VDS with
    @return: VDS with chosen annotations and filters
    """
    vds = hl.vds.read_vds(vds_path, n_partitions=n_partitions)

    if hard_filtered_samples_to_remove_ht:
        vds = hl.vds.filter_samples(vds, hard_filtered_samples_to_remove_ht, keep=False)

    if release_only:
        assert meta_ht is not None
        meta_ht = meta_ht.filter(meta_ht.release)
        vds = hl.vds.filter_samples(vds, meta_ht)

    if split:
        vmt = vds.variant_data
        vmt = vmt.annotate_rows(
            n_unsplit_alleles=hl.len(vmt.alleles),
            mixed_site=(hl.len(vmt.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(vmt.alleles[0], a), vmt.alleles[1:])
            & hl.any(lambda a: hl.is_snp(vmt.alleles[0], a), vmt.alleles[1:]),
        )
        vmt = hl.experimental.sparse_split_multi(vmt, filter_changed_loci=True)
        vds = hl.vds.VariantDataset(vds.reference_data, vmt)

    return vds
