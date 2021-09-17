"""
Dataproc cluster utils
"""
from typing import Optional, List
from analysis_runner import dataproc
from hailtop.batch.job import Job
from hailtop.batch import Batch
from joint_calling import utils


clusters_by_name = dict()


def get_cluster(
    b: Batch,
    name,
    num_workers: int,
    long: bool = False,
    is_test: bool = False,
    phantomjs: bool = False,
    preempt: bool = True,
    depends_on: Optional[List[Job]] = None,
):
    """
    Get or create a Dataproc cluster by name
    """
    if name not in clusters_by_name:
        max_age = '1h' if is_test else '12h'
        if long:
            max_age = '2h' if is_test else '24h'

        clusters_by_name[name] = dataproc.setup_dataproc(
            b,
            max_age=max_age,
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=num_workers if preempt else 0,
            num_workers=2 if preempt else num_workers,
            cluster_name=name,
            depends_on=depends_on,
            init=['gs://cpg-reference/hail_dataproc/install_phantomjs.sh']
            if phantomjs
            else [],
        )
    return clusters_by_name[name]
