import collections
import logging

from larcoh.pipeline import Pipeline
from larcoh.targets import Cohort
from larcoh.tasks.combiner import Combiner
from larcoh.types import SequencingType

logger = logging.getLogger(__file__)


def subset_to_test(
    cohort: Cohort,
    dataset: str | None = None,
):
    """
    Function that subsets the cohort for test runs. By default, takes 30 samples.
    Also take optional dataset name for custom subset code.
    """
    if dataset == 'prophecy':
        # Test dataset: collecting all batch1 samples, plus 4 samples from each of other batches
        samples_by_batch: dict[str, list] = collections.defaultdict(list)
        for sample in cohort.get_samples():
            batch = sample.sequencing_meta_by_type[SequencingType.GENOME]['batch']
            print(
                sample, sample.sequencing_meta_by_type[SequencingType.GENOME]['batch']
            )
            if batch != '1' and len(samples_by_batch[batch]) > 4:
                sample.active = False
            samples_by_batch[batch].append(sample)


def main():
    pipeline = Pipeline()

    pipeline.add_task(Combiner(pipeline))

    # relatedness_prefix = pipeline.analysis_prefix / 'relatedness'
    #
    # sample_qc_job, hard_filter_ht_path, meta_ht_path = add_sample_qc_jobs(
    #     b=b,
    #     mt_path=raw_combined_mt_path,
    #     samples_tsv_path=samples_tsv_path,
    #     sample_qc_bucket=join(analysis_bucket, 'sample_qc'),
    #     ancestry_bucket=join(analysis_bucket, 'ancestry'),
    #     tmp_bucket=tmp_bucket,
    #     analysis_bucket=analysis_bucket,
    #     relatedness_bucket=relatedness_bucket,
    #     release_related=release_related,
    #     web_bucket=web_bucket,
    #     filter_cutoffs_path=filter_cutoffs_path,
    #     overwrite=overwrite,
    #     scatter_count=scatter_count,
    #     combiner_job=combiner_job,
    #     billing_project=billing_project,
    #     is_test=output_namespace in ['test', 'tmp'],
    #     highmem_workers=highmem_workers,
    # )
    #
    # var_qc_job = add_variant_qc_jobs(
    #     b=b,
    #     work_bucket=join(analysis_bucket, 'variant_qc'),
    #     web_bucket=join(web_bucket, 'variant_qc'),
    #     raw_combined_mt_path=raw_combined_mt_path,
    #     hard_filter_ht_path=hard_filter_ht_path,
    #     meta_ht_path=meta_ht_path,
    #     out_filtered_combined_mt_path=filtered_combined_mt_path,
    #     sample_count=len(samples_df),
    #     ped_file=ped_fpath,
    #     overwrite=overwrite,
    #     vqsr_params_d=utils.get_filter_cutoffs(filter_cutoffs_path)['vqsr'],
    #     scatter_count=scatter_count,
    #     is_test=output_namespace in ['test', 'tmp'],
    #     depends_on=[combiner_job, sample_qc_job],
    #     highmem_workers=highmem_workers,
    # )
    #
    # job_name = 'Remove ref blocks'
    # if not utils.can_reuse(filtered_combined_noref_mt_path, overwrite):
    #     noref_mt_j = add_job(
    #         b,
    #         f'{utils.SCRIPTS_DIR}/make_noref_mt.py '
    #         f'--overwrite '
    #         f'--mt {filtered_combined_mt_path} '
    #         f'--out-mt {filtered_combined_noref_mt_path}',
    #         job_name=job_name,
    #         num_workers=scatter_count,
    #         depends_on=[var_qc_job],
    #     )
    # else:
    #     noref_mt_j = b.new_job(f'{job_name} [reuse]')
    #
    # if subset_projects:
    #     diff_projects = set(subset_projects) - set(input_projects)
    #     if diff_projects:
    #         raise click.BadParameter(
    #             f'--subset-project values should be a subset of --input-project '
    #             f'values. The following projects are not in input projects: '
    #             f'{diff_projects} '
    #         )
    #     subset_projects = list(set(subset_projects))
    #     subset_mt_path = f'{release_bucket}/mt/{output_version}-{"-".join(sorted(subset_projects))}.mt'
    #     job_name = f'Making subset MT for {", ".join(subset_projects)}'
    #     if overwrite or not utils.file_exists(subset_mt_path):
    #         add_job(
    #             b,
    #             f'{utils.SCRIPTS_DIR}/make_subset_mt.py '
    #             f'--mt {filtered_combined_noref_mt_path} ' +
    #             (''.join(f'--subset-project {p} ' for p in subset_projects)) +
    #             f'--out-mt {subset_mt_path}',
    #             job_name=job_name,
    #             is_test=output_namespace in ['test', 'tmp'],
    #             num_workers=scatter_count,
    #             depends_on=[noref_mt_j],
    #         )
    #     else:
    #         b.new_job(f'{job_name} [reuse]')
    #
    #     subset_vcf_path = f'{release_bucket}/vcf/{output_version}-{"-".join(subset_projects)}.vcf.bgz'
    #     job_name = f'Convert subset MT to VCF for {", ".join(subset_projects)}'
    #     if overwrite or not utils.file_exists(subset_vcf_path):
    #         add_job(
    #             b,
    #             f'{utils.SCRIPTS_DIR}/final_mt_to_vcf.py '
    #             f'--mt {subset_mt_path} ' +
    #             f'--out-vcf {subset_vcf_path}',
    #             job_name=job_name,
    #             is_test=output_namespace in ['test', 'tmp'],
    #             num_workers=scatter_count,
    #             preemptible=False,
    #             depends_on=[noref_mt_j],
    #         )
    #     else:
    #         b.new_job(f'{job_name} [reuse]')

    pipeline.run()


if __name__ == '__main__':
    main()
