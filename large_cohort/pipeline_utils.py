"""
Utility functions for the pipeline orchestration code. Depend on cpg-pipes 
and python 3.10, hence not importable in Dataproc scripts.
"""

import logging
from functools import lru_cache

import math
from analysis_runner import dataproc
from cpg_pipes.hb.batch import RegisteringBatch, setup_batch
from cpg_pipes.utils import timestamp, slugify
from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path
from hailtop.batch import Batch
from hailtop.batch.job import Job

logger = logging.getLogger(__file__)


DATAPROC_PACKAGES = [
    'cpg-utils',
    'coloredlogs',
    'click',
    'cpg-gnomad==0.6.3',
    'google',
    'slackclient',
    'fsspec',
    'sklearn',
    'gcloud',
    'selenium',
]


@lru_cache
def get_batch(cohort) -> RegisteringBatch:
    """
    Lazy initialise Batch object.
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


# GCP quota on number of cores
MAX_PRIMARY_WORKERS = 50


def dataproc_job(
    batch: Batch,
    script_name: str,
    params: dict,
    preemptible: bool = True,
    phantomjs: bool = True,
    num_workers: int | None = None,
    depends_on: list[Job | None] = None,
    autoscaling_policy: str | None = None,
    long: bool = False,
    worker_boot_disk_size: int | None = None,
    secondary_worker_boot_disk_size: int | None = None,
) -> Job:
    """
    Submit script as a dataproc job.
    """
    param_line = ' '.join(
        f'--{k.replace("_", "-")} {v}' for k, v in params.items() if v is not None
    )
    if num_workers is None:
        num_workers = get_config()['workflow']['scatter_count']

    max_age = '24h' if not long else '48h'

    depends_on = depends_on or []
    depends_on = [j for j in depends_on if j is not None]

    if preemptible:
        num_secondary_workers = num_workers
        # number of primary workers has to be 5-10% of the number of secondary workers:
        # see Tim's comment in https://discuss.hail.is/t/all-nodes-are-unhealthy/1764
        # using 8% here:
        num_primary_workers = int(math.ceil(num_secondary_workers * 0.08))
    else:
        num_secondary_workers = 0
        num_primary_workers = num_workers

    use_highmem_workers = get_config()['workflow'].get('highmem_workers')

    # 2 is the minimal number of primary workers for dataproc cluster:
    num_primary_workers = max(num_primary_workers, 2)

    # limiting the number of primary workers to avoid hitting GCP quota:
    num_primary_workers = min(num_primary_workers, MAX_PRIMARY_WORKERS)

    return dataproc.hail_dataproc_job(
        batch,
        script=f'scripts/{script_name} {param_line}',
        job_name=script_name,
        max_age=max_age,
        packages=DATAPROC_PACKAGES,
        num_secondary_workers=num_secondary_workers,
        num_workers=num_primary_workers,
        autoscaling_policy=autoscaling_policy,
        depends_on=depends_on,
        init=['gs://cpg-reference/hail_dataproc/install_phantomjs.sh']
        if phantomjs
        else [],
        worker_machine_type='n1-highmem-8' if use_highmem_workers else 'n1-standard-8',
        worker_boot_disk_size=worker_boot_disk_size,
        secondary_worker_boot_disk_size=secondary_worker_boot_disk_size,
        pyfiles=['large_cohort'],
    )
