"""
Functions to find the pipeline inputs and communicate with the SM server
"""

import logging
import sys
from dataclasses import dataclass
from os.path import join
from typing import List, Dict, Optional, Set, Collection

import pandas as pd
from hailtop.batch import Batch
from hailtop.batch.job import Job

from sample_metadata.apis import (
    AnalysisApi,
    SequenceApi,
    SampleApi,
)
from sample_metadata.models import (
    AnalysisStatus,
    AnalysisType,
    AnalysisQueryModel,
)
from sample_metadata import exceptions

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
    timestamp_completed: str
    meta: Dict


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
        timestamp_completed=data.get('timestamp_completed', None),
        meta=data.get('meta', {}),
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
    project: str,
    analysis_type: str,
    analysis_status: str = 'completed',
    meta: Optional[Dict] = None,
) -> Dict[str, Analysis]:
    """
    Query the DB to find the last completed analysis for the type and samples,
    one Analysis object per sample. Assumes the analysis is defined for a single
    sample (e.g. cram, gvcf)
    """
    analysis_per_sid: Dict[str, Analysis] = dict()

    datas = aapi.query_analyses(
        AnalysisQueryModel(
            projects=[project],
            type=AnalysisType(analysis_type),
            status=AnalysisStatus(analysis_status),
            sample_ids=sample_ids,
            meta=meta or {},
        )
    )

    for data in datas:
        a = _parse_analysis(data)
        if not a:
            continue
        assert a.type == analysis_type, data
        assert a.status == analysis_status, data
        assert len(a.sample_ids) == 1, data
        sid = list(a.sample_ids)[0]
        analysis_per_sid[sid] = a

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


default_entry = {
    's': None,
    'external_id': None,
    'stack': None,
    'project': None,
    'source': '-',
    'continental_pop': '-',
    'subpop': '-',
    'topostproc_gvcf': '-',
    'gvcf': '-',
    'cram': '-',
    'crai': '-',
    'realign_cram': '-',
    'realign_crai': '-',
    'batch': '-',
    'resequencing_label': '-',
    'primary_study': '-',
    'operation': 'add',
    'flowcell_lane': '-',
    'library_id': '-',
    'platform': '-',
    'centre': '-',
    'r_contamination': None,
    'r_chimera': None,
    'r_duplication': None,
    'median_insert_size': None,
    'median_coverage': None,
    'r_30x': None,
    'r_aligned_in_pairs': None,
    'fam_id': '-',
    'mat_id': '-',
    'pat_id': '-',
    'sex': '-',
    'sex_karyotype': '-',
    'age': None,
}


def find_inputs_from_db(
    input_projects: List[str],
    skip_samples: Optional[Collection[str]] = None,
    check_existence: bool = True,
    source_tag: Optional[str] = None,
    assume_gvcfs_are_ready: bool = False,
) -> pd.DataFrame:
    """
    Determine input samples and pull input files and metadata from the CPG
    sample-metadata server database.

    To specify a subset of GVCFs to search, set `source_tag`, which should
    match the meta.source value analysis entries.
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

        meta = {}
        if source_tag:
            meta['source'] = source_tag

        qc_analysis_per_sid = find_analyses_by_sid(
            sample_ids=[s['id'] for s in samples],
            analysis_type='qc',
            project=proj,
            meta=dict(**meta),
        )

        reblocked_gvcf_by_sid: Dict[str, str] = {}
        staging_gvcf_by_sid: Dict[str, str] = {}
        if assume_gvcfs_are_ready:
            if proj.endswith('-test'):
                stack = proj.replace('-test', '')
                namespace = 'test'
            else:
                stack = proj
                namespace = 'main'
            proj_bucket = f'gs://cpg-{stack}-{namespace}/gvcf'
            if source_tag:
                proj_bucket = join(proj_bucket, source_tag)
            reblocked_gvcf_by_sid = {
                s['id']: join(proj_bucket, f'{s["id"]}.g.vcf.gz')
                for s in samples
            }
        else:
            reblocked_gvcf_analysis_per_sid = find_analyses_by_sid(
                sample_ids=[s['id'] for s in samples],
                analysis_type='gvcf',
                project=proj,
                meta=dict(staging=False, **meta),
            )
            staging_gvcf_analysis_per_sid = dict()
            if source_tag is not None:
                staging_gvcf_analysis_per_sid = find_analyses_by_sid(
                    sample_ids=[s['id'] for s in samples],
                    project=proj,
                    analysis_type='gvcf',
                    meta=dict(staging=True, **meta),
                )

            sids_without_gvcf = []
            for s in samples:
                reblocked_gvcf_analysis = reblocked_gvcf_analysis_per_sid.get(s['id'])
                staging_gvcf_analysis = staging_gvcf_analysis_per_sid.get(s['id'])
                if not reblocked_gvcf_analysis and not staging_gvcf_analysis:
                    sids_without_gvcf.append(s['id'] + '/' + s['external_id'])
                    continue
    
                if reblocked_gvcf_analysis:
                    if not reblocked_gvcf_analysis.output:
                        logger.error(
                            f'"output" is not defined for the latest reblocked gvcf '
                            f' analysis, skipping sample {s["id"]}'
                        )
                        sids_without_gvcf.append(s['id'] + '/' + s['external_id'])
                        continue
                    reblocked_gvcf_by_sid[s['id']] = reblocked_gvcf_analysis.output
    
                if staging_gvcf_analysis:
                    if not staging_gvcf_analysis.output:
                        logger.error(
                            f'"output" is not defined for the latest staging gvcf analysis, '
                            f'skipping sample {s["id"]}'
                        )
                        sids_without_gvcf.append(s['id'] + '/' + s['external_id'])
                        continue
                    staging_gvcf_by_sid[s['id']] = staging_gvcf_analysis.output
    
            if sids_without_gvcf:
                logger.warning(
                    f'No gvcf found for {len(sids_without_gvcf)}/{len(samples)} samples: '
                    f'{", ".join(sids_without_gvcf)}'
                )
            samples = [
                s
                for s in samples
                if s['id'] in reblocked_gvcf_by_sid or s['id'] in staging_gvcf_by_sid
            ]
            if not samples:
                logger.info(f'No samples to process, skipping project {proj}')
                continue

        logger.info('Checking sequencing info for samples')
        seq_meta_by_sid = dict()
        try:
            seq_infos: List[Dict] = seqapi.get_sequences_by_sample_ids(
                request_body=[s['id'] for s in samples]
            )
        except exceptions.ApiException:
            logger.warning(f'Sequencing data was not found for some samples')
        else:
            sids_without_meta = []
            for si in seq_infos:
                if 'meta' not in si:
                    sids_without_meta.append(si['sample_id'])
                else:
                    seq_meta_by_sid[si['sample_id']] = si['meta']
            if sids_without_meta:
                logger.warning(
                    f'Found {len(sids_without_meta)} samples without "meta" in '
                    f'sequencing info: {", ".join(sids_without_meta)}'
                )
            logger.info(f'Found {len(seq_meta_by_sid)} sequences')

        logger.info(
            'Checking GVCFs for samples and collecting pipeline input DataFrame'
        )
        for s in samples:
            sample_id = s['id']
            external_id = s['external_id']
            seq_meta = seq_meta_by_sid.get(sample_id, {})
            reblocked_gvcf_path = reblocked_gvcf_by_sid.get(s['id'])
            staging_gvcf_path = staging_gvcf_by_sid.get(s['id'])

            qc_analysis = qc_analysis_per_sid.get(sample_id)
            if qc_analysis:
                qc_metrics = qc_analysis.meta.get('metrics', {})
            else:
                qc_metrics = {} 

            def _check_gvcf(gvcf_path):
                if not gvcf_path.endswith('.g.vcf.gz'):
                    logger.warning(
                        f'GVCF analysis for sample ID {sample_id} "output" field '
                        f'is not a GVCF'
                    )
                    return False
                if check_existence and not utils.file_exists(gvcf_path):
                    logger.warning(
                        f'GVCF analysis for sample ID {sample_id} "output" file '
                        f'does not exist: {gvcf_path}'
                    )
                    return False
                if check_existence and not utils.file_exists(gvcf_path + '.tbi'):
                    logger.warning(
                        f'GVCF analysis for sample ID {sample_id} "output" field '
                        f'does not have a corresponding tbi index: {gvcf_path}.tbi'
                    )
                    return False
                return True

            if not _check_gvcf(reblocked_gvcf_path or staging_gvcf_path):
                continue
                
            stack = proj.replace('-test', '')
            project = s['meta'].get('project', proj.replace('-test', ''))

            entry = default_entry.copy()
            entry.update(
                {
                    's': sample_id,
                    'external_id': external_id,
                    'fam_id': external_id,
                    'stack': stack,
                    'project': project,
                    'source': source_tag or '-',
                    'gvcf': reblocked_gvcf_path or '-',
                    'topostproc_gvcf': staging_gvcf_path or '-',
                    'batch': seq_meta.get('batch', '-'),
                    'flowcell_lane': seq_meta.get(
                        'sample.flowcell_lane', seq_meta.get('flowcell_lane', '-')
                    ),
                    'library_id': seq_meta.get(
                        'sample.library_id', seq_meta.get('library_id', '-')
                    ),
                    'platform': seq_meta.get(
                        'sample.platform', seq_meta.get('platform', '-')
                    ),
                    'centre': seq_meta.get(
                        'sample.centre', seq_meta.get('centre', '-')
                    ),
                    'primary_study': seq_meta.get('Primary study', '-'),
                    # QC metrics:
                    'r_contamination': qc_metrics.get('freemix'),
                    'r_chimera': qc_metrics.get('pct_chimeras'),
                    'r_duplication': qc_metrics.get('percent_duplication'),
                    'median_insert_size': qc_metrics.get('median_insert_size'),
                    'median_coverage': qc_metrics.get('median_coverage'),
                    'r_30x': qc_metrics.get('pct_30x'),
                    'r_aligned_in_pairs': qc_metrics.get('pct_reads_aligned_in_pairs'),
                    'continental_pop': s['meta'].get('continental_pop', '-'),
                    'subcontinental_pop': s['meta'].get('subcontinental_pop') or s['meta'].get('subpop', '-'),
                }
            )
            inputs.append(entry)

    if not inputs:
        logger.error('Noе found any projects with samples good for processing')
        sys.exit(1)

    df = pd.DataFrame(inputs).set_index('s', drop=False)
    return df


def add_validation_samples(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add NA12878 GVCFs and syndip BAM into the dataframe.
    """
    if 'syndip' not in df.s:
        entry = default_entry.copy()
        entry.update(
            {
                's': 'syndip',
                'external_id': 'syndip',
                'fam_id': 'syndip',
                'project': 'syndip',
                'cram': 'gs://cpg-reference/validation/syndip/raw/CHM1_CHM13_2.bam',
                'crai': 'gs://cpg-reference/validation/syndip/raw/CHM1_CHM13_2.bam.bai',
            }
        )
        # Can only append a dict if ignore_index=True. So then need to set index back.
        df = df.append(entry, ignore_index=True).set_index('s', drop=False)

    giab_samples = ['NA12878', 'NA12891', 'NA12892']
    for sn in giab_samples:
        if sn not in df.s:
            cram = f'gs://cpg-reference/validation/giab/cram/{sn}.cram'
            entry = default_entry.copy()
            entry.update(
                {
                    's': sn,
                    'external_id': sn,
                    'fam_id': 'CEPH',
                    'sex': '1' if sn == 'NA12891' else '2',
                    'sex_karyotype': 'XY' if sn == 'NA12891' else 'XX',
                    'pat_id': 'NA12891' if sn == 'NA12878' else '-',
                    'mat_id': 'NA12892' if sn == 'NA12878' else '-',
                    'project': 'giab',
                    'cram': cram,
                    'crai': cram + '.crai',
                }
            )
            # Can only append a dict if ignore_index=True. So then need to set index back.
            df = df.append(entry, ignore_index=True).set_index('s', drop=False)
    return df


@dataclass
class AlignmentInput:
    """
    Sort of a union type for possible alignment inputs
    """

    bam_or_cram_path: Optional[str] = None
    index_path: Optional[str] = None
    fqs1: Optional[List[str]] = None
    fqs2: Optional[List[str]] = None
