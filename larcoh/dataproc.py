"""
Dataproc cluster utils
"""
import math
from typing import Optional, List

from analysis_runner import dataproc

from hailtop.batch.job import Job
from hailtop.batch import Batch

from larcoh import utils


DataprocCluster = dataproc.DataprocCluster


# GCP quota on number of cores
MAX_PRIMARY_WORKERS = 50


def add_job(
    b: Batch,
    script: str,
    job_name: str,
    num_workers: int,
    long: bool = False,
    highmem: bool = False,
    is_test: bool = False,
    phantomjs: bool = False,
    preemptible: bool = True,
    depends_on: Optional[List[Optional[Job]]] = None,
    autoscaling_policy: Optional[str] = None,
    max_age: Optional[str] = None,
    pyfiles: Optional[List[str]] = None,
    worker_boot_disk_size: Optional[int] = None,
    secondary_worker_boot_disk_size: Optional[int] = None,
):
    """
    Wrapper around submitting a dataproc cluster job
    """
    if not max_age:
        max_age = '1h' if is_test else '24h'
        if long:
            max_age = '3h' if is_test else '48h'

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

    # 2 is the minimal number of primary workers for dataproc cluster:
    num_primary_workers = max(num_primary_workers, 2)

    # limiting the number of primary workers to avoid hitting GCP quota:
    num_primary_workers = min(num_primary_workers, MAX_PRIMARY_WORKERS)
    
    return dataproc.hail_dataproc_job(
        b,
        script=script,
        job_name=job_name,
        max_age=max_age,
        packages=utils.DATAPROC_PACKAGES,
        num_secondary_workers=num_secondary_workers,
        num_workers=num_primary_workers,
        autoscaling_policy=autoscaling_policy,
        depends_on=depends_on,
        init=['gs://cpg-reference/hail_dataproc/install_phantomjs.sh']
        if phantomjs
        else [],
        worker_machine_type='n1-highmem-8' 
        if highmem 
        else 'n1-standard-8',
        pyfiles=pyfiles,
        worker_boot_disk_size=worker_boot_disk_size,
        secondary_worker_boot_disk_size=secondary_worker_boot_disk_size,
    )


def get_cluster(
    b: Batch,
    name,
    num_workers: int,
    long: bool = False,
    highmem: bool = False,
    is_test: bool = False,
    phantomjs: bool = False,
    preemptible: bool = True,
    depends_on: Optional[List[Optional[Job]]] = None,
    autoscaling_policy: Optional[str] = None,
    max_age: Optional[str] = None,
) -> dataproc.DataprocCluster:
    """
    Get or create a Dataproc cluster by name
    """
    if not max_age:
        max_age = '1h' if is_test else '24h'
        if long:
            max_age = '3h' if is_test else '48h'

    depends_on = depends_on or []
    depends_on = [j for j in depends_on if j is not None]

    return dataproc.setup_dataproc(
        b,
        max_age=max_age,
        packages=utils.DATAPROC_PACKAGES,
        num_secondary_workers=num_workers if preemptible else 0,
        num_workers=2 if preemptible else num_workers,
        cluster_name=name,
        depends_on=depends_on,
        init=['gs://cpg-reference/hail_dataproc/install_phantomjs.sh']
        if phantomjs
        else [],
        autoscaling_policy=autoscaling_policy,
        worker_machine_type='n1-highmem-8' 
        if highmem 
        else 'n1-standard-8',
    )
