# """
# Sample-level variant QC and sex imputation
# """
#
# from cpg_pipes.hb.batch import RegisteringBatch
# from cpg_pipes.targets import Cohort
# from cpg_utils import Path
# from hailtop.batch.job import Job
#
# from larcoh.pipeline_utils import dataproc_job
# from larcoh.query_utils import can_reuse
#
#
# def sample_qc(
#     batch: RegisteringBatch,
#     cohort: Cohort,
#     vds_path: Path,
#     out_ht_path: Path,
#     depends_on: list[Job] | None = None,
# ) -> Job | None:
