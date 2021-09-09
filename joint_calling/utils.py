"""Utility functions for the joint_calling module"""

import os
import subprocess
import tempfile
import logging
import sys
import time
import hashlib
from os.path import isdir, isfile, exists, join
from typing import Any, Callable, Dict, Optional
import yaml
import hail as hl
import click
from google.cloud import storage
from joint_calling import _version, get_package_path


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


DEFAULT_REF = 'GRCh38'

REF_BUCKET = 'gs://cpg-reference/hg38/v1'

DATAPROC_PACKAGES = [
    'joint-calling',
    'click',
    'cpg-gnomad',
    'google',
    'slackclient',
    'fsspec',
    'sklearn',
    'gcloud',
]

DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver'

AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
GATK_VERSION = '4.2.0.0'
GATK_IMAGE = f'{AR_REPO}/gatk:{GATK_VERSION}'
# GnarlyGenotyper is in Beta and crashes with NullPointerException when using the
# official GATK docker, that's why we're using a separate image for it:
GNARLY_IMAGE = f'{AR_REPO}/gnarly_genotyper:hail_ukbb_300K'
BCFTOOLS_IMAGE = f'{AR_REPO}/bcftools:1.10.2--h4f4756c_2'
SM_IMAGE = f'{AR_REPO}/sm-api:2.0.3'

NOALT_REGIONS = join(REF_BUCKET, 'noalt.bed')
TEL_AND_CENT_HT_PATH = join(
    REF_BUCKET, 'gnomad/telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht'
)
LCR_INTERVALS_HT_PATH = join(REF_BUCKET, 'gnomad/lcr_intervals/LCRFromHengHg38.ht')
SEG_DUP_INTERVALS_HT_PATH = join(
    REF_BUCKET, 'gnomad/seg_dup_intervals/GRCh38_segdups.ht'
)
CLINVAR_HT_PATH = join(REF_BUCKET, 'gnomad/clinvar/clinvar_20190923.ht')

GNOMAD_HT_PATH = (
    'gs://gcp-public-data--gnomad/release/3.1/ht/genomes/gnomad.genomes.v3.1.sites.ht/'
)
GNOMAD_HGDP_1KG_MT_PATH = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)
GNOMAD_HGDP_1KG_TEST_MT_PATH = join(
    REF_BUCKET, 'mt/gnomad.genomes.v3.1.hgdp_1kg_subset_dense_test.mt'
)


NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS = 50
NUMBER_OF_GENOMICS_DB_INTERVALS = 50


# Sample QC output file names
POP_HT_NAME = 'pop.ht'
REGRESSED_METRICS_HT_NAME = 'regressed_metrics.ht'
RELATED_SAMPLES_TO_DROP_HT_NAME = 'related_samples_to_drop.ht'

PCA_EIGENVALUES_HT_NAME = 'pca_eigenvalues.ht'
PCA_SCORES_HT_NAME = 'pca_scores.ht'
PCA_LOADINGS_HT_NAME = 'pca_loadings.ht'


def init_hail(name: str, local_tmp_dir: str = None):
    """
    Initialize Hail and set up the directory for logs
    :param name: name to prefix the log file
    :param local_tmp_dir: local directory to write Hail logs
    :return:
    """
    if not local_tmp_dir:
        local_tmp_dir = tempfile.mkdtemp()

    timestamp = time.strftime('%Y%m%d-%H%M')
    hl_log = os.path.join(
        safe_mkdir(os.path.join(local_tmp_dir, 'log')), f'{name}-{timestamp}.log'
    )
    hl.init(default_reference=DEFAULT_REF, log=hl_log)
    logger.info(f'Running joint-calling version {_version.__version__}')
    return local_tmp_dir


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
          in which case it will check for the existence of a
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


def can_reuse(fpath: Optional[str], overwrite: bool, silent=False) -> bool:
    """
    Checks if the file `fpath` exists and we are not overwriting
    """
    if not fpath:
        return False
    if not file_exists(fpath):
        return False

    if overwrite:
        if not silent:
            logger.info(f'File {fpath} exists and will be overwritten')
        return False
    else:
        if not silent:
            logger.info(f'Reusing existing {fpath}. Use --overwrite to overwrite')
        return True


def gs_cache_file(fpath: str, local_tmp_dir: str) -> str:
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


def get_mt(
    mt_path: str,
    split: bool = False,
    hard_filtered_samples_to_remove_ht: hl.Table = None,
    meta_ht: hl.Table = None,
    add_meta: bool = False,
    release_only: bool = False,
    passing_sites_only: bool = False,
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
        release (can only be used if metadata is present)
    :param passing_sites_only: whether to filter the MT to only variants with
        nothing in the filter field (e.g. passing soft filters)
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

    elif release_only:
        assert meta_ht is not None
        mt = mt.filter_cols(meta_ht[mt.col_key].release)

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


def get_vqsr_filters_path(
    work_bucket: str,
    model_id: str,
    split: bool = True,
    finalized: bool = False,
) -> str:
    """
    Gets the specified VQSR filtering annotation resource.
    :param work_bucket: bucket
    :param model_id: VQSR filtering model id
    :param split: Split or multi-allelic version of the filtering file
    :param finalized: Whether to return the raw VQSR table or the finalized VQSR table representing determined cutoffs
    :return: VQSR filtering annotation file path
    """
    return join(
        work_bucket,
        f'filtering/{model_id}'
        f'{".finalized" if finalized else ""}'
        f'{".split" if split else ""}'
        f'.ht',
    )


def get_filter_cutoffs(
    provided_filter_cutoffs_path: Optional[str] = None,
) -> Dict:
    """
    :provided_filter_cutoffs_path: optional, a path to a YAML file with cutoffs.
    Can sit on a bucket. If not provided, a default one from the package will be used.
    gets the a default one within the package
    :return: a Dict with cutoffs
    """
    if provided_filter_cutoffs_path:
        assert file_exists(provided_filter_cutoffs_path), provided_filter_cutoffs_path
        path = provided_filter_cutoffs_path
    else:
        path = join(get_package_path(), 'filter_cutoffs.yaml')

    if path.startswith('gs://'):
        contents = subprocess.check_output(['gsutil', 'cat', path])
        filter_cutoffs_d = yaml.safe_load(contents)
    else:
        with open(path) as f:
            filter_cutoffs_d = yaml.safe_load(f)

    return filter_cutoffs_d
