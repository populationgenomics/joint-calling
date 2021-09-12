"""
Submit a hail query script that pre-subsets the HGDP matrix table
"""

import os
import hailtop.batch as hb
from analysis_runner import dataproc


POP = 'nfe'


service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'),
    bucket=os.getenv('HAIL_BUCKET'),
)
batch = hb.Batch(name='subset-gnomad-hgdp-1kg', backend=service_backend)

make_full_hq_j = dataproc.hail_dataproc_job(
    batch,
    'subset_hgdp.py',
    max_age='12h',
    packages=['selenium'],
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name='subset-gnomad-hgdp-1kg',
    num_secondary_workers=50,
)

make_test_hq_j = dataproc.hail_dataproc_job(
    batch,
    'subset_hgdp.py --test',
    max_age='12h',
    packages=['selenium'],
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name='subset-gnomad-hgdp-1kg',
    num_secondary_workers=50,
    depends_on=[make_full_hq_j],
)

make_pop_hq_j = dataproc.hail_dataproc_job(
    batch,
    'subset_hgdp.py --pop',
    max_age='12h',
    packages=['selenium'],
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name='subset-gnomad-hgdp-1kg',
    num_secondary_workers=50,
    depends_on=[make_full_hq_j],
)

dataproc.hail_dataproc_job(
    batch,
    f'subset_hgdp.py --test --pop {POP}',
    max_age='12h',
    packages=['selenium'],
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    job_name='subset-gnomad-hgdp-1kg',
    num_secondary_workers=50,
    depends_on=[make_test_hq_j],
)

batch.run()
