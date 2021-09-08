"""
Functions to find the pipeline inputs and communicate with the SM server
"""

import logging
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from os.path import join, basename
from typing import List, Dict, Optional, Set, Collection

import pandas as pd
from hailtop.batch import Batch
from hailtop.batch.job import Job
from sample_metadata import (
    AnalysisApi,
    SequenceApi,
    SampleApi,
    exceptions,
)
from joint_calling import utils

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


sapi = SampleApi()
aapi = AnalysisApi()
seqapi = SequenceApi()


@dataclass
class Analysis:
    """
    Represents the analysis DB entry
    """

    id: str
    type: str
    status: str
    sample_ids: Set[str]
    output: Optional[str]


def _parse_analysis(data: Dict) -> Optional[Analysis]:
    if not data:
        return None
    if 'id' not in data:
        logger.error(f'Analysis data doesn\'t have id: {data}')
        return None
    if 'type' not in data:
        logger.error(f'Analysis data doesn\'t have type: {data}')
        return None
    if 'status' not in data:
        logger.error(f'Analysis data doesn\'t have status: {data}')
        return None
    a = Analysis(
        id=data['id'],
        type=data['type'],
        status=data['status'],
        sample_ids=set(data.get('sample_ids', [])),
        output=data.get('output', None),
    )
    return a


def find_joint_calling_analysis(
    analysis_project: str,
    sample_ids: Collection[str],
) -> Optional[Analysis]:
    """
    Query the DB to find the last completed joint-calling analysis for the samples
    """
    data = aapi.get_latest_complete_analysis_for_type(
        project=analysis_project,
        analysis_type='joint-calling',
    )
    a = _parse_analysis(data)
    if not a:
        return None
    assert a.type == 'joint-calling', data
    assert a.status == 'completed', data
    if a.sample_ids != set(sample_ids):
        return None
    return a


def find_analyses_by_sid(
    sample_ids: Collection[str],
    analysis_project: str,
    analysis_type: str,
) -> Dict[str, Analysis]:
    """
    Query the DB to find the last completed analysis for the type and samples,
    one Analysis object per sample. Assumes the analysis is defined for a single
    sample (e.g. cram, gvcf)
    """
    analysis_per_sid: Dict[str, Analysis] = dict()
    datas = aapi.get_latest_analysis_for_samples_and_type(
        project=analysis_project,
        analysis_type=analysis_type,
        request_body=sample_ids,
    )
    for data in datas:
        a = _parse_analysis(data)
        if not a:
            continue
        assert a.type == analysis_type, data
        assert a.status == 'completed', data
        assert len(a.sample_ids) == 1, data
        analysis_per_sid[list(a.sample_ids)[0]] = a
    return analysis_per_sid


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
    j.image(utils.SM_IMAGE)
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


def replace_paths_to_test(s: Dict) -> Optional[Dict]:
    """
    Replace paths of all files in -main namespace to -test namespsace,
    and return None if files in -test are not found.
    :param s:
    :return:
    """

    def fix(fpath):
        fpath = fpath.replace('-main-upload/', '-test-upload/')
        if not utils.file_exists(fpath):
            return None
        return fpath

    try:
        reads_type = s['meta']['reads_type']
        if reads_type in ('bam', 'cram'):
            fpath = s['meta']['reads'][0]['location']
            fpath = fix(fpath)
            if not fpath:
                return None
            s['meta']['reads'][0]['location'] = fpath

            fpath = s['meta']['reads'][0]['secondaryFiles'][0]['location']
            fpath = fix(fpath)
            if not fpath:
                return None
            s['meta']['reads'][0]['secondaryFiles'][0]['location'] = fpath

        elif reads_type == 'fastq':
            for li in range(len(s['meta']['reads'])):
                for rj in range(len(s['meta']['reads'][li])):
                    fpath = s['meta']['reads'][li][rj]['location']
                    fpath = fix(fpath)
                    if not fpath:
                        return None
                    s['meta']['reads'][li][rj]['location'] = fpath

        logger.info(f'Found test sample {s["id"]}')
        return s
    except Exception:  # pylint: disable=broad-except
        return None


def find_inputs_from_db(
    input_projects: List[str],
    is_test: bool = False,
    skip_samples: Optional[Collection[str]] = None,
) -> pd.DataFrame:
    """
    Determine inputs from SM DB
    """
    inputs = []

    for proj in input_projects:
        logger.info(f'Processing project {proj}')
        samples = sapi.get_samples(
            body_get_samples_by_criteria_api_v1_sample_post={
                'project_ids': [proj],
                'active': True,
            }
        )
        logger.info(f'Found {len(samples)} samples')
        if not samples:
            logger.info(f'No samples to process, skipping project {proj}')
            continue

        if skip_samples:
            logger.info('Checking which samples need to skip')
            not_skipped_sids = []
            for s in samples:
                if skip_samples and s['id'] in skip_samples:
                    logger.info(f'Skiping sample: {s["id"]}')
                    continue
                not_skipped_sids.append(s['id'])
            logger.info(f'Excluding skipped samples: {len(not_skipped_sids)}')
            samples = [s for s in samples if s in not_skipped_sids]
            if not samples:
                logger.info(f'No samples to process, skipping project {proj}')
                continue

        logger.info('Checking GVCF analyses for samples')
        gvcf_analysis_per_sid = find_analyses_by_sid(
            sample_ids=[s['id'] for s in samples],
            analysis_type='gvcf',
            analysis_project=proj,
        )
        gvcf_by_sid = dict()
        sids_without_gvcf = []
        for s in samples:
            a = gvcf_analysis_per_sid.get(s['id'])
            if not a:
                sids_without_gvcf.append(s['id'])
                continue
            gvcf_path = a.output
            if not gvcf_path:
                logger.error(
                    f'"output" is not defined for the latest gvcf analysis, '
                    f'skipping sample {s["id"]}'
                )
                sids_without_gvcf.append(s['id'])
                continue
            gvcf_by_sid[s['id']] = gvcf_path

        if sids_without_gvcf:
            logger.warning(
                f'No gvcf found for {len(sids_without_gvcf)}/{len(samples)} samples: '
                f'{", ".join(sids_without_gvcf)}'
            )
            samples = [s for s in samples if s['id'] in gvcf_by_sid]
            if not samples:
                logger.info(f'No samples to process, skipping project {proj}')
                continue

        logger.info('Checking sequencing info for samples')
        try:
            seq_infos: List[Dict] = seqapi.get_sequences_by_sample_ids(
                request_body=[s['id'] for s in samples]
            )
        except exceptions.ApiException:
            logger.critical(f'Not for all samples sequencing data was found')
            raise

        seq_meta_by_sid: Dict = dict()
        sids_without_meta = []
        for si in seq_infos:
            if 'meta' not in si:
                sids_without_meta.append(si['sample_id'])
            else:
                seq_meta_by_sid[si['sample_id']] = si['meta']
        if sids_without_meta:
            logger.error(
                f'Found {len(sids_without_meta)} samples without "meta" in '
                f'sequencing info: {", ".join(sids_without_meta)}'
            )
        samples = [s for s in samples if s['id'] in seq_meta_by_sid]
        if not samples:
            logger.info(f'No samples to process, skipping project {proj}')
            continue

        logger.info(
            'Checking GVCFs for samples and collecting pipeline input data frame'
        )
        for s in samples:
            sample_id = s['id']
            external_id = s['external_id']
            seq_meta = seq_meta_by_sid[sample_id]
            gvcf_path = gvcf_by_sid[s['id']]

            # TODO: reenable once we support row data and crams
            # if is_test:
            #     s = replace_paths_to_test(s)
            # if s:
            #     samples_by_project[proj].append(s)

            if is_test:
                if '/batch1/' not in gvcf_path:
                    continue
                gvcf_path = gvcf_path.replace(
                    f'gs://cpg-{proj}-main',
                    f'gs://cpg-{proj}-test',
                )
                gvcf_path = gvcf_path.replace(s['id'], s['external_id'])
                if not utils.file_exists(gvcf_path):
                    continue
                logger.info(f'Using {gvcf_path} for a test run')

            if not gvcf_path.endswith('.g.vcf.gz'):
                logger.warning(
                    f'GVCF analysis for sample ID {sample_id} "output" field '
                    f'is not a GVCF'
                )
                continue
            if not utils.file_exists(gvcf_path):
                logger.warning(
                    f'GVCF analysis for sample ID {sample_id} "output" file '
                    f'does not exist: {gvcf_path}'
                )
                continue
            if not utils.file_exists(gvcf_path + '.tbi'):
                logger.warning(
                    f'GVCF analysis for sample ID {sample_id} "output" field '
                    f'does not have a corresponding tbi index: {gvcf_path}.tbi'
                )
                continue
            sample_information = {
                's': sample_id,
                'external_id': external_id,
                'population': None,
                'gvcf': gvcf_path,
                'batch': seq_meta.get('batch'),
                'r_contamination': seq_meta.get('raw_data.FREEMIX'),
                'r_chimera': seq_meta.get('raw_data.PCT_CHIMERAS'),
                'r_duplication': seq_meta.get('raw_data.PERCENT_DUPLICATION'),
                'median_insert_size': seq_meta.get('raw_data.MEDIAN_INSERT_SIZE'),
                'flowcell_lane': seq_meta.get('sample.flowcell_lane'),
                'library_id': seq_meta.get('sample.library_id'),
                'platform': seq_meta.get('sample.platform'),
                'centre': seq_meta.get('sample.centre'),
                'operation': 'add',
            }
            inputs.append(sample_information)

    if not inputs:
        logger.error('No found any projects with samples good for processing')
        sys.exit(1)

    df = pd.DataFrame(inputs).set_index('s', drop=False)
    return df


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
            single_df['population'] = None
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
                population=None,
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
