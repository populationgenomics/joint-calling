import importlib.metadata
from typing import Callable

import coloredlogs
import logging
import collections
from functools import lru_cache

from cpg_utils.config import get_config
from cpg_utils import Path, to_path
from cpg_utils.hail_batch import dataset_path

from hailtop.batch.job import Job

from larcoh.dataproc import add_job
from larcoh.hb.batch import setup_batch, RegisteringBatch
from larcoh.utils import can_reuse, timestamp, slugify, exists
from larcoh.filetypes import SequencingType
from larcoh.inputs.metamist.inputs import CpgInputProvider
from larcoh.inputs.metamist.metamist import Metamist
from larcoh.targets import Cohort


coloredlogs.install(
    level='DEBUG', fmt='%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
)


logger = logging.getLogger(__file__)


def get_package_name() -> str:
    """
    Get name of the package.
    """
    return __name__.split('.', 1)[0]


def get_package_path() -> Path:
    """
    Get local install path of the package.
    """
    return to_path(__file__).parent.absolute()


def get_version() -> str:
    """
    Get package version.
    """
    return importlib.metadata.version(get_package_name())


output_version = get_config()['workflow']['output_version']

analysis_prefix = to_path(
    dataset_path(f'joint-calling/{output_version}', category='analysis')
)

tmp_prefix = to_path(dataset_path(f'joint-calling/{output_version}', category='tmp'))


@lru_cache
def cohort() -> Cohort:
    """
    Read inputs from the database.
    """
    c = Cohort(
        analysis_dataset_name=get_config()['workflow']['dataset'],
        sequencing_type=SequencingType.parse(
            get_config()['workflow']['sequencing_type']
        ),
    )

    CpgInputProvider(db=Metamist(c.analysis_dataset.name)).populate_cohort(
        cohort=c,
        dataset_names=get_config()['workflow'].get('datasets'),
        skip_samples=get_config()['workflow'].get('skip_samples'),
        only_samples=get_config()['workflow'].get('only_samples'),
        skip_datasets=get_config()['workflow'].get('skip_datasets'),
    )

    if get_config()['workflow'].get('access_level') == 'test':
        # Function that subsets the cohort for test runs
        original_sample_cnt = len(c.get_samples())
        original_dataset_cnt = len(c.get_datasets())
        # Test dataset: collecting all batch1 samples, plus 4 samples from each of other batches
        samples_by_batch: dict[str, list] = collections.defaultdict(list)
        for sample in c.get_samples():
            batch_id = sample.sequencing_meta_by_type[SequencingType.GENOME]['batch']
            if batch != '1' and len(samples_by_batch[batch_id]) > 4:
                sample.active = False
                continue
            samples_by_batch[batch_id].append(sample)
        logger.info(
            f'After subset: {len(c.get_samples())}/{original_sample_cnt} samples, '
            f'{len(samples_by_batch)} batches, '
            f'{len(c.get_datasets())}/{original_dataset_cnt} datasets'
        )
        logger.info(
            f'Using {len(c.get_samples())} samples '
            f'from {len(c.get_datasets())} datasets'
        )

    return c


@lru_cache
def batch() -> RegisteringBatch:
    """
    Initialise Batch object
    """
    run_id = get_config()['workflow'].get('run_id', timestamp())
    analysis_dataset = get_config()['workflow']['dataset']
    name = get_config()['workflow'].get('name')
    description = get_config()['workflow'].get('description')
    name = name or description or analysis_dataset
    name = slugify(name)
    description = description or name
    description += f': run_id={run_id}'

    description += f' [{cohort().sequencing_type.value}]'
    if ds_set := set(d.name for d in cohort().get_datasets()):
        description += ' ' + ', '.join(sorted(ds_set))
    return setup_batch(description=description)


def dataproc_job(
    script_name: str,
    params: dict,
    preemptible: bool = True,
    phantomjs: bool = True,
    num_workers: int | None = None,
    depends_on: list[Job] = None,
    autoscaling_policy: str | None = None,
    long: bool = False,
) -> Job:
    param_line = ' '.join(
        f'--{k.replace("_", "-")} {v}' for k, v in params.items() if v is not None
    )
    if num_workers is None:
        num_workers = get_config()['workflow']['scatter_count']

    return add_job(
        batch(),
        f'scripts/{script_name} {param_line}',
        job_name=script_name,
        num_workers=num_workers,
        highmem=get_config()['workflow'].get('highmem_workers'),
        preemptible=preemptible,
        phantomjs=phantomjs,
        depends_on=depends_on or [],
        autoscaling_policy=autoscaling_policy,
        long=long,
    )
