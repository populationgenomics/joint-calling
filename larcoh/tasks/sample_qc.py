from cpg_utils import Path
from hailtop.batch.job import Job

from larcoh import batch, utils, can_reuse, dataproc_job, cohort, tmp_prefix


def sample_qc(
    mt_path: Path,
    sex_ht_path: Path,
    hard_filtered_samples_ht_path: Path,
    hail_sample_qc_ht_path: Path,
    depends_on: list[Job] | None = None,
) -> list[Job]:
    if can_reuse(
        [
            hard_filtered_samples_ht_path,
            sex_ht_path,
            hail_sample_qc_ht_path,
        ]
    ):
        sample_qc_hardfilter_jobs = [batch().new_job(f'Sample QC [reuse]')]
    else:
        sample_qc_hardfilter_jobs = dataproc_job(
            'sample_qc_hard_filters.py',
            params=dict(
                mt=mt_path,
                cohort_tsv=cohort().to_tsv(),
                out_hard_filtered_samples_ht_path=hard_filtered_samples_ht_path,
                out_sex_ht_path=sex_ht_path,
                out_hail_sample_qc_ht_path=hail_sample_qc_ht_path,
                tmp_prefix=tmp_prefix / 'sample_qc',
            ),
            depends_on=depends_on,
        )
    return sample_qc_hardfilter_jobs
