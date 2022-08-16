"""
Utility functions for the pipeline orchestration code. Depend on cpg-pipes 
and python 3.10, hence not importable in Dataproc scripts.
"""

import logging
from functools import lru_cache

from cpg_pipes.hb.batch import RegisteringBatch, setup_batch
from cpg_pipes.utils import timestamp, slugify
from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path
from hailtop.batch.job import Job

from larcoh import dataproc

logger = logging.getLogger(__file__)


DATAPROC_PACKAGES = [
    'larcoh',
    'click',
    'cpg-gnomad',
    'google',
    'slackclient',
    'fsspec',
    'sklearn',
    'gcloud',
    'selenium',
]


output_version = get_config()['workflow']['output_version']
vds_version = get_config()['workflow'].get('vds_version', output_version)

_suffix = f'joint-calling/{output_version}'
analysis_prefix = to_path(dataset_path(_suffix, category='analysis'))
tmp_prefix = to_path(dataset_path(_suffix, category='tmp'))


@lru_cache
def get_batch(cohort) -> RegisteringBatch:
    """
    Initialise Batch object
    """
    run_id = get_config()['workflow'].get('run_id', timestamp())
    analysis_dataset = get_config()['workflow']['dataset']
    name = get_config()['workflow'].get('name')
    description = get_config()['workflow'].get('description')
    sequencing_type = get_config()['workflow']['sequencing_type']
    name = name or description or analysis_dataset
    name = slugify(name)
    description = description or name
    description += f': run_id={run_id} [{sequencing_type}]'
    if ds_set := set(d.name for d in cohort.get_datasets()):
        description += ' ' + ', '.join(sorted(ds_set))
    return setup_batch(description=description)


def dataproc_job(
    batch: RegisteringBatch,
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

    return dataproc.add_job(
        batch,
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
