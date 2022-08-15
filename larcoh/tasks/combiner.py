import logging

from cpg_utils import Path
from cpg_utils.config import get_config
from hailtop.batch.job import Job

from larcoh import (
    dataproc_job,
    tmp_prefix,
)
from larcoh.utils import exists, can_reuse

logger = logging.getLogger(__file__)


def queue_combiner(batch, cohort, out_vds_path: Path) -> Job | None:
    """
    Add VCF combiner jobs, produce VDC.
    """
    if can_reuse(out_vds_path):
        return None

    # Combiner takes advantage of autoscaling cluster policies
    # to reduce costs for the work that uses only the driver machine:
    # https://hail.is/docs/0.2/experimental/vcf_combiner.html#pain-points
    # To add a 50-worker policy for a project "prophecy-339301":
    # ```
    # gcloud dataproc autoscaling-policies import vcf-combiner-50 \
    # --source=combiner-autoscaling-policy-50.yaml --region=australia-southeast1 \
    # --project prophecy-339301
    # ```
    scatter_count = get_config()['workflow'].get('scatter_count', 50)
    if scatter_count > 100:
        autoscaling_workers = '200'
    elif scatter_count > 50:
        autoscaling_workers = '100'
    else:
        autoscaling_workers = '50'

    for i, sample in enumerate(cohort.get_samples()):
        if not sample.gvcf:
            if get_config()['workflow'].get('skip_samples_with_missing_input', False):
                logger.warning(f'Skipping {sample} which is missing GVCF')
                sample.active = False
                continue
            else:
                raise ValueError(
                    f'Sample {sample} is missing GVCF. '
                    f'Use workflow/skip_samples = [] or '
                    f'workflow/skip_samples_with_missing_input '
                    f'to control behaviour'
                )

        gvcf_path = sample.gvcf.path
        if get_config()['workflow'].get('check_inputs', True):
            if not exists(gvcf_path):
                if get_config()['workflow'].get(
                    'skip_samples_with_missing_input', False
                ):
                    logger.warning(
                        f'Skipping {sample} that is missing GVCF {gvcf_path}'
                    )
                    sample.active = False
                else:
                    raise ValueError(
                        f'Sample {sample} is missing GVCF. '
                        f'Use workflow/skip_samples = [] or '
                        f'workflow/skip_samples_with_missing_input '
                        f'to control behaviour'
                    )

    logger.info(
        f'Combining {len(cohort.get_samples())} samples: '
        f'{", ".join(cohort.get_sample_ids())}'
    )

    branch_factor = get_config().get('combiner', {}).get('branch_factor')
    batch_size = get_config().get('combiner', {}).get('batch_size')

    job = dataproc_job(
        batch,
        'combine_gvcfs.py',
        params=dict(
            cohort_tsv=cohort.to_tsv(),
            out_vds=out_vds_path,
            tmp_prefix=tmp_prefix / 'combiner',
            branch_factor=branch_factor,
            batch_size=batch_size,
        ),
        num_workers=0,
        autoscaling_policy=f'vcf-combiner-{autoscaling_workers}',
        long=True,
    )
    return job
