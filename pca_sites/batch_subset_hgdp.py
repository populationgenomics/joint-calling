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

make_pop_hq_j = dataproc.hail_dataproc_job(
    batch,
    f'subset_hgdp.py --pop {POP}',
    max_age='12h',
    packages=['selenium'],
    job_name=f'subset-gnomad-hgdp-1kg {POP}',
    num_secondary_workers=50,
)


batch.run()
