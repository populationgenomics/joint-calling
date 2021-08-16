"""Utility functions for the joint_calling module"""

import os
import subprocess
import tempfile
import logging
import sys
import time
import hashlib
from dataclasses import dataclass
from os.path import isdir, isfile, exists, join, basename
from typing import Any, Callable, List, Dict, Optional, Tuple
import shutil
import yaml
import pandas as pd
import hail as hl
import click
from google.cloud import storage
from hailtop.batch import Batch
from hailtop.batch.job import Job
from sample_metadata import (
    AnalysisApi,
    AnalysisType,
    AnalysisStatus,
)
from joint_calling import _version, get_package_path

logger = logging.getLogger('joint-calling')
logger.setLevel('INFO')


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

TEL_AND_CENT_HT_PATH = join(
    REF_BUCKET, 'gnomad/telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht'
)
LCR_INTERVALS_HT_PATH = join(REF_BUCKET, 'gnomad/lcr_intervals/LCRFromHengHg38.ht')
SEG_DUP_INTERVALS_HT_PATH = join(
    REF_BUCKET, 'gnomad/seg_dup_intervals/GRCh38_segdups.ht'
)
CLINVAR_HT_PATH = join(REF_BUCKET, 'gnomad/clinvar/clinvar_20190923.ht')


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


@dataclass
class Analysis:
    """
    Represents the analysis DB entry
    """

    id: str
    type: AnalysisType
    status: AnalysisStatus
    output: Optional[str]
    sample_ids: List[str]

    @staticmethod
    def from_db(**kwargs):
        """
        Convert from db keys, mainly converting id to id_
        """
        analysis_type = kwargs.pop('type', None)
        status = kwargs.pop('status', None)
        sample_ids = kwargs['sample_ids']
        output = kwargs.pop('output', [])
        return Analysis(
            id=kwargs.pop('id'),
            type=AnalysisType(analysis_type),
            status=AnalysisStatus(status),
            sample_ids=list(set(sorted(sample_ids))),
            output=output,
        )


def get_latest_complete_analysis(
    analysis_project: str,
) -> Dict[Tuple[str, Tuple], Analysis]:
    """
    Returns a dictionary that maps a tuple (analysis type, sample ids) to the
    lastest complete analysis record (represented by a AnalysisModel object)
    """
    aapi = AnalysisApi()
    latest_by_type_and_sids = dict()
    for a_type in ['cram', 'gvcf', 'joint-calling']:
        for a_data in aapi.get_latest_complete_analyses_by_type(
            project=analysis_project,
            analysis_type=a_type,
        ):
            a: Analysis = Analysis.from_db(**a_data)
            latest_by_type_and_sids[(a_type, tuple(a.sample_ids))] = a
    return latest_by_type_and_sids


def make_sm_in_progress_job(
    b: Batch, analyais_type: str, analysis_project: str, analysis_id: str
) -> Job:
    """
    Creates a job that updates the sample metadata server entry analysis status
    to in-progress
    """
    return make_sm_update_status_job(
        b, analyais_type, 'in-progress', analysis_project, analysis_id
    )


def make_sm_completed_job(
    b: Batch, analyais_type: str, sm_db_name: str, analysis_id: str
) -> Job:
    """
    Creates a job that updates the sample metadata server entry analysis status
    to completed
    """
    return make_sm_update_status_job(
        b, analyais_type, 'completed', sm_db_name, analysis_id
    )


def make_sm_update_status_job(
    b: Batch, analysis_type: str, status: str, sm_db_name: str, analysis_id: str
) -> Job:
    """
    Creates a job that updates the sample metadata server entry analysis status.
    """
    assert status in ['in-progress', 'failed', 'completed', 'queued']
    j = b.new_job(f'SM: update {analysis_type} to {status}')
    j.image(SM_IMAGE)
    j.command(
        f"""
set -o pipefail
set -ex

export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS
export SM_DEV_DB_PROJECT={sm_db_name}
export SM_ENVIRONMENT=PRODUCTION

cat <<EOT >> update.py
from sample_metadata.api import AnalysisApi
from sample_metadata import AnalysisUpdateModel
aapi = AnalysisApi()
aapi.update_analysis_status(
    analysis_id='{analysis_id}',
    analysis_update_model=AnalysisUpdateModel(status='{status}'),
)
EOT
python update.py
    """
    )
    return j


def find_inputs(
    input_buckets: List[str],
    input_metadata_buckets: Optional[List[str]] = None,
) -> pd.DataFrame:  # pylint disable=too-many-branches
    """
    Read the inputs assuming a standard CPG storage structure.
    :param input_buckets: buckets to find GVCFs
    :param input_metadata_buckets: buckets to find CSV metadata files
    :return: a dataframe with the following structure:
        s (key)
        population
        gvcf
        r_contamination
        r_chimera
        r_duplication
        median_insert_size
    """
    gvcf_paths: List[str] = []
    for ib in input_buckets:
        cmd = f'gsutil ls \'{ib}/*.g.vcf.gz\''
        gvcf_paths.extend(
            line.strip()
            for line in subprocess.check_output(cmd, shell=True).decode().split()
        )

    local_tmp_dir = tempfile.mkdtemp()

    if input_metadata_buckets:
        qc_csvs: List[str] = []
        for ib in input_metadata_buckets:
            cmd = f'gsutil ls \'{ib}/*.csv\''
            qc_csvs.extend(
                line.strip()
                for line in subprocess.check_output(cmd, shell=True).decode().split()
            )

        df: pd.DataFrame = None
        # sample.id,sample.sample_name,sample.flowcell_lane,sample.library_id,sample.platform,sample.centre,sample.reference_genome,raw_data.FREEMIX,raw_data.PlinkSex,raw_data.PCT_CHIMERAS,raw_data.PERCENT_DUPLICATION,raw_data.MEDIAN_INSERT_SIZE,raw_data.MEDIAN_COVERAGE
        # 613,TOB1529,ILLUMINA,HVTVGDSXY.1-2-3-4,LP9000039-NTP_H04,KCCG,hg38,0.0098939700,F(-1),0.023731,0.151555,412.0,31.0
        # 609,TOB1653,ILLUMINA,HVTVGDSXY.1-2-3-4,LP9000039-NTP_F03,KCCG,hg38,0.0060100100,F(-1),0.024802,0.165634,452.0,33.0
        # 604,TOB1764,ILLUMINA,HVTV7DSXY.1-2-3-4,LP9000037-NTP_B02,KCCG,hg38,0.0078874400,F(-1),0.01684,0.116911,413.0,43.0
        # 633,TOB1532,ILLUMINA,HVTVGDSXY.1-2-3-4,LP9000039-NTP_C05,KCCG,hg38,0.0121946000,F(-1),0.024425,0.151094,453.0,37.0
        columns = {
            'sample.sample_name': 's',
            'raw_data.FREEMIX': 'r_contamination',
            'raw_data.PCT_CHIMERAS': 'r_chimera',
            'raw_data.PERCENT_DUPLICATION': 'r_duplication',
            'raw_data.MEDIAN_INSERT_SIZE': 'median_insert_size',
        }
        for qc_csv in qc_csvs:
            local_qc_csv_path = join(local_tmp_dir, basename(qc_csv))
            subprocess.run(
                f'gsutil cp {qc_csv} {local_qc_csv_path}', check=False, shell=True
            )
            single_df = pd.read_csv(local_qc_csv_path)
            single_df = single_df.rename(columns=columns)[columns.values()]
            single_df['population'] = 'EUR'
            single_df['gvcf'] = ''
            single_df = single_df.set_index('s', drop=False)
            df = (
                single_df
                if df is None
                else (pd.concat([df, single_df], ignore_index=False).drop_duplicates())
            )
        sample_names = list(df['s'])
    else:
        sample_names = [basename(gp).replace('.g.vcf.gz', '') for gp in gvcf_paths]
        df = pd.DataFrame(
            data=dict(
                s=sample_names,
                population='EUR',
                gvcf=gvcf_paths,
                r_contamination=pd.NA,
                r_chimera=pd.NA,
                r_duplication=pd.NA,
                median_insert_size=pd.NA,
            )
        ).set_index('s', drop=False)

    shutil.rmtree(local_tmp_dir)

    # Checking 1-to-1 match of sample names to GVCFs
    for sn in sample_names:
        matching_gvcfs = [gp for gp in gvcf_paths if sn in gp]
        if len(matching_gvcfs) > 1:
            logging.warning(
                f'Multiple GVCFs found for the sample {sn}:' f'{matching_gvcfs}'
            )
        elif len(matching_gvcfs) == 0:
            logging.warning(f'No GVCFs found for the sample {sn}')

    # Checking 1-to-1 match of GVCFs to sample names, and filling a dict
    for gp in gvcf_paths:
        matching_sn = [sn for sn in sample_names if sn in gp]
        if len(matching_sn) > 1:
            logging.warning(
                f'Multiple samples found for the GVCF {gp}:' f'{matching_sn}'
            )
        elif len(matching_sn) == 0:
            logging.warning(f'No samples found for the GVCF {gp}')
        else:
            df.loc[matching_sn[0], ['gvcf']] = gp
    df = df[df.gvcf.notnull()]
    return df


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
        filter_cutoffs_d = yaml.load(contents)
    else:
        with open(path) as f:
            filter_cutoffs_d = yaml.load(f)

    return filter_cutoffs_d
