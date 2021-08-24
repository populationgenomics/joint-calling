"""
Functions to find the pipeline inputs and communicate with the SM server
"""

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from os.path import join, basename
from typing import List, Dict, Optional, Tuple

import pandas as pd
from hailtop.batch import Batch
from hailtop.batch.job import Job
from sample_metadata import (
    AnalysisApi,
    AnalysisType,
    AnalysisStatus,
    SequenceApi,
    SampleApi,
)
from joint_calling import utils

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


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


def find_inputs_from_db(
    input_projects: List[str], analysis_project: str, is_test: bool = False
) -> pd.DataFrame:
    """
    Determine inputs from SM DB
    """
    aapi = AnalysisApi()
    seqapi = SequenceApi()
    sapi = SampleApi()

    # Get samples in latest analysis
    latest_analysis = aapi.get_latest_complete_analyses_by_type(
        analysis_type='joint-calling', project=analysis_project
    )
    latest_analysis_sample_ids = [
        analysis['sample_ids'] for analysis in latest_analysis
    ]
    sample_ids_in_latest_analysis = [
        sample_id for sublist in latest_analysis_sample_ids for sample_id in sublist
    ]

    active_samples = sapi.get_samples(
        body_get_samples_by_criteria_api_v1_sample_post={
            'project_ids': input_projects,
            'active': True,
        }
    )
    active_samples_by_id = {s['id']: s for s in active_samples}
    active_sample_ids = set(active_samples_by_id.keys())
    new_sample_ids = list(set(active_sample_ids) - set(sample_ids_in_latest_analysis))
    deleted_sample_ids = list(
        set(sample_ids_in_latest_analysis) - set(active_sample_ids)
    )

    inputs = []

    # Get all sequence metadata for the list of processed samples
    sequences_data = seqapi.get_sequences_by_sample_ids(request_body=new_sample_ids)

    new_sample_gvcf_analyses = aapi.get_latest_gvcfs_for_samples(new_sample_ids)

    for new_gvcf_analysis in new_sample_gvcf_analyses:
        sample_id = new_gvcf_analysis.get('sample_ids')[0]
        current_seq_data = next(
            seq_data
            for seq_data in sequences_data
            if seq_data['sample_id'] == sample_id
        )
        gvcf_path = new_gvcf_analysis.get('output')
        external_id = active_samples_by_id[sample_id]['external_id']
        if is_test:
            if '/batch1/' not in gvcf_path:
                continue
            for proj in input_projects:
                gvcf_path = gvcf_path.replace(
                    f'gs://cpg-{proj}-main',
                    f'gs://cpg-{proj}-test',
                )
                if not utils.file_exists(gvcf_path):
                    continue
            gvcf_path = gvcf_path.replace(sample_id, external_id)
            logger.info(f'Using {gvcf_path} for a test run')

        if not gvcf_path:
            logger.warning(
                f'GVCF analysis for sample ID {sample_id} does not have '
                f'an "output" field'
            )
            continue
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
            'population': 'EUR',
            'gvcf': gvcf_path,
            'r_contamination': current_seq_data.get('meta').get('raw_data.FREEMIX'),
            'r_chimera': current_seq_data.get('meta').get('raw_data.PCT_CHIMERAS'),
            'r_duplication': current_seq_data.get('meta').get(
                'raw_data.PERCENT_DUPLICATION'
            ),
            'median_insert_size': current_seq_data.get('meta').get(
                'raw_data.MEDIAN_INSERT_SIZE'
            ),
            'operation': 'add',
        }

        inputs.append(sample_information)

    for sample_id in deleted_sample_ids:
        sample_information = {
            's': sample_id,
            'external_id': None,
            'population': None,
            'gvcf': None,
            'r_contamination': None,
            'r_chimera': None,
            'r_duplication': None,
            'median_insert_size': None,
            'operation': 'delete',
        }
        inputs.append(sample_information)

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
