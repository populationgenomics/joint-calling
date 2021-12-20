"""
Dataproc cluster utils
"""
from typing import Optional, List
from analysis_runner import dataproc
from hailtop.batch.job import Job
from hailtop.batch import Batch
from joint_calling import utils


DataprocCluster = dataproc.DataprocCluster


# GCP quota on number of cores
MAX_NONPREEMPTIBLE_WORKERS = 50


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
):
    if not max_age:
        max_age = '1h' if is_test else '24h'
        if long:
            max_age = '3h' if is_test else '48h'

    depends_on = depends_on or []
    depends_on = [j for j in depends_on if j is not None]
    
    return dataproc.hail_dataproc_job(
        b,
        script=script,
        job_name=job_name,
        max_age=max_age,
        packages=utils.DATAPROC_PACKAGES,
        num_secondary_workers=num_workers if preemptible else 0,
        num_workers=min(2 if preemptible else num_workers, MAX_NONPREEMPTIBLE_WORKERS),
        autoscaling_policy=autoscaling_policy,
        depends_on=depends_on,
        init=['gs://cpg-reference/hail_dataproc/install_phantomjs.sh']
        if phantomjs
        else [],
        worker_machine_type='n1-highmem-8' 
        if highmem 
        else 'n1-standard-8',
        pyfiles=pyfiles,
    )


def get_cluster(
    b: Batch,
    name,
    num_workers: int,
    long: bool = False,
    highmem_workers: bool = False,
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
        if highmem_workers 
        else 'n1-standard-8',
    )
