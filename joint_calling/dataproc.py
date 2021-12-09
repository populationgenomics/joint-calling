"""
Dataproc cluster utils
"""
from typing import Optional, List
from analysis_runner import dataproc
from hailtop.batch.job import Job
from hailtop.batch import Batch
from joint_calling import utils


DataprocCluster = dataproc.DataprocCluster


def get_cluster(
    b: Batch,
    name,
    num_workers: int,
    long: bool = False,
    is_test: bool = False,
    phantomjs: bool = False,
    preemptible: bool = True,
    depends_on: Optional[List[Optional[Job]]] = None,
) -> dataproc.DataprocCluster:
    """
    Get or create a Dataproc cluster by name
    """
    max_age = '1h' if is_test else '8h'
    if long:
        max_age = '3h' if is_test else '24h'
    
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
    )
