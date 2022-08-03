"""
Helpers to communicate with the sample-metadata database.
"""

import logging
import pprint
from dataclasses import dataclass
from enum import Enum

from cpg_utils import to_path
from sample_metadata.apis import (
    SampleApi,
    SequenceApi,
    AnalysisApi,
    ParticipantApi,
    FamilyApi,
)

from larcoh import utils
from larcoh.filetypes import (
    FastqPair,
    CramPath,
    AlignmentInput,
    SequencingType,
    FastqPairs,
)

logger = logging.getLogger(__file__)


class SmdbError(Exception):
    """
    Raised for problems interacting with sample-metadata database.
    """


class AnalysisType(Enum):
    """
    Corresponds to SMDB Analysis types:
    https://github.com/populationgenomics/sample-metadata/blob/dev/models/enums
    /analysis.py#L4-L11

    Re-defined in a separate module to decouple from the main sample-metadata module,
    so decorators can use `@stage(analysis_type=AnalysisType.QC)` without importing
    the sample-metadata package.
    """

    QC = 'qc'
    JOINT_CALLING = 'joint-calling'
    GVCF = 'gvcf'
    CRAM = 'cram'
    CUSTOM = 'custom'
    ES_INDEX = 'es-index'

    @staticmethod
    def parse(val: str) -> 'AnalysisType':
        """
        Parse str and create a AnalysisStatus object
        """
        d = {v.value: v for v in AnalysisType}
        if val not in d:
            raise SmdbError(
                f'Unrecognised analysis type {val}. Available: {list(d.keys())}'
            )
        return d[val.lower()]


class Metamist:
    """
    Communication with the SampleMetadata database.
    """

    def __init__(self, project_name: str | None = None):
        """
        @param project_name: default SMDB project name.
        """
        self.sapi = SampleApi()
        self.aapi = AnalysisApi()
        self.seqapi = SequenceApi()
        self.seqapi = SequenceApi()
        self.papi = ParticipantApi()
        self.fapi = FamilyApi()
        self.project_name = project_name

    def get_sample_entries(
        self,
        project_name: str | None = None,
        active: bool = True,
    ) -> list[dict]:
        """
        Get samples in the project as a list of dictionaries.
        """
        project_name = project_name or self.project_name

        logger.debug(f'Finding samples for dataset {project_name}...')
        body = {
            'project_ids': [project_name],
            'active': active,
        }
        sample_entries = self.sapi.get_samples(body_get_samples=body)
        logger.info(
            f'Finding samples for project {project_name}: '
            f'found {len(sample_entries)}'
        )
        return sample_entries

    def get_ped_entries(self, project_name: str | None = None) -> list[dict[str, str]]:
        """
        Retrieve PED lines for a specified SM project, with external participant IDs.
        """
        project_name = project_name or self.project_name

        families = self.fapi.get_families(project_name)
        family_ids = [family['id'] for family in families]
        ped_entries = self.fapi.get_pedigree(
            internal_family_ids=family_ids,
            response_type='json',
            project=project_name,
            replace_with_participant_external_ids=True,
        )

        return ped_entries


@dataclass
class MmSequence:
    """
    Sample-metadata DB "Sequence" entry.

    See sample-metadata for more details:
    https://github.com/populationgenomics/sample-metadata
    """

    id: str
    sample_id: str
    meta: dict
    sequencing_type: SequencingType
    alignment_input: AlignmentInput | None = None

    @staticmethod
    def parse(data: dict, check_existence: bool) -> 'MmSequence':
        """
        Parse dictionary to create a SmSequence object.
        """
        req_keys = ['id', 'sample_id', 'meta']
        if any(k not in data for k in req_keys):
            for key in req_keys:
                if key not in data:
                    logger.error(f'"Sequence" data does not have {key}: {data}')
            raise ValueError(f'Cannot parse SMDB Sequence {data}')

        sample_id = data['sample_id']
        st = SequencingType.parse(data['type'])
        assert st, data
        sm_seq = MmSequence(
            id=data['id'],
            sample_id=sample_id,
            meta=data['meta'],
            sequencing_type=st,
        )
        if data['meta'].get('reads'):
            if alignment_input := MmSequence._parse_reads(
                sample_id=sample_id,
                meta=data['meta'],
                check_existence=check_existence,
            ):
                sm_seq.alignment_input = alignment_input
        else:
            logger.warning(
                f'{sample_id} sequence: no meta/reads found with FASTQ information'
            )
        return sm_seq

    @staticmethod
    def _parse_reads(  # pylint: disable=too-many-return-statements
        sample_id: str,
        meta: dict,
        check_existence: bool,
    ) -> AlignmentInput | None:
        """
        Parse a AlignmentInput object from the meta dictionary.

        @param check_existence: check if fastq/crams exist on buckets.
        Default value is pulled from self.smdb and can be overridden.
        """
        reads_data = meta.get('reads')
        reads_type = meta.get('reads_type')
        reference_assembly = meta.get('reference_assembly', {}).get('location')

        if not reads_data:
            logger.error(f'{sample_id}: no "meta/reads" field in meta')
            return None
        if not reads_type:
            logger.error(f'{sample_id}: no "meta/reads_type" field in meta')
            return None
        supported_types = ('fastq', 'bam', 'cram')
        if reads_type not in supported_types:
            logger.error(
                f'{sample_id}: ERROR: "reads_type" is expected to be one of '
                f'{supported_types}'
            )
            return None

        if reads_type in ('bam', 'cram'):
            if len(reads_data) > 1:
                logger.error(f'{sample_id}: supporting only single bam/cram input')
                return None

            bam_path = reads_data[0]['location']
            if not (bam_path.endswith('.cram') or bam_path.endswith('.bam')):
                logger.error(
                    f'{sample_id}: ERROR: expected the file to have an extension '
                    f'.cram or .bam, got: {bam_path}'
                )
                return None
            if check_existence and not utils.exists(bam_path):
                logger.error(
                    f'{sample_id}: ERROR: index file does not exist: {bam_path}'
                )
                return None

            # Index:
            index_path = None
            if reads_data[0].get('secondaryFiles'):
                index_path = reads_data[0]['secondaryFiles'][0]['location']
                if (
                    bam_path.endswith('.cram')
                    and not index_path.endswith('.crai')
                    or bam_path.endswith('.bai')
                    and not index_path.endswith('.bai')
                ):
                    logger.error(
                        f'{sample_id}: ERROR: expected the index file to have an extension '
                        f'.crai or .bai, got: {index_path}'
                    )
                if check_existence and not utils.exists(index_path):
                    logger.error(
                        f'{sample_id}: ERROR: index file does not exist: {index_path}'
                    )
                    return None

            return CramPath(
                bam_path, index_path=index_path, reference_assembly=reference_assembly
            )

        else:
            fastq_pairs = FastqPairs()
            for lane_pair in reads_data:
                if len(lane_pair) != 2:
                    raise ValueError(
                        f'Sequence data for sample {sample_id} is incorrectly '
                        f'formatted. Expecting 2 entries per lane (R1 and R2 fastqs), '
                        f'but got {len(lane_pair)}. '
                        f'Read data: {pprint.pformat(reads_data)}'
                    )
                if check_existence and not utils.exists(lane_pair[0]['location']):
                    logger.error(
                        f'{sample_id}: ERROR: read 1 file does not exist: '
                        f'{lane_pair[0]["location"]}'
                    )
                    return None
                if check_existence and not utils.exists(lane_pair[1]['location']):
                    logger.error(
                        f'{sample_id}: ERROR: read 2 file does not exist: '
                        f'{lane_pair[1]["location"]}'
                    )
                    return None

                fastq_pairs.append(
                    FastqPair(
                        to_path(lane_pair[0]['location']),
                        to_path(lane_pair[1]['location']),
                    )
                )

            return fastq_pairs
