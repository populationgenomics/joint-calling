import collections
import logging

from cpg_pipes.inputs import get_cohort
from cpg_pipes.targets import Cohort
from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path

from larcoh import (
    analysis_prefix,
    get_batch,
    vds_version,
)
from larcoh.tasks import combiner
from larcoh.tasks.sample_qc import sample_qc

logger = logging.getLogger(__file__)


def _test_subset(c: Cohort) -> Cohort:
    """
    Subset cohort the cohort for test runs
    """
    original_sample_cnt = len(c.get_samples())
    original_dataset_cnt = len(c.get_datasets())
    # Test dataset: collecting all batch1 samples, plus 4 samples from each of other batches
    samples_by_batch: dict[str, list] = collections.defaultdict(list)
    for sample in c.get_samples():
        batch_id = sample.seq_by_type['genome'].meta['batch']
        if get_batch != '1' and len(samples_by_batch[batch_id]) > 4:
            sample.active = False
            continue
        samples_by_batch[batch_id].append(sample)
    logger.info(
        f'After subset: {len(c.get_samples())}/{original_sample_cnt} samples, '
        f'{len(samples_by_batch)} batches, '
        f'{len(c.get_datasets())}/{original_dataset_cnt} datasets'
    )
    logger.info(
        f'Using {len(c.get_samples())} samples '
        f'from {len(c.get_datasets())} datasets'
    )
    return c


def main():
    cohort = get_cohort()
    if get_config()['workflow'].get('access_level') == 'test':
        cohort = _test_subset(cohort)

    batch = get_batch(cohort)

    # COMBINER
    vds_path = to_path(dataset_path(f'{vds_version}.vds'))
    combiner_job = combiner.queue_combiner(batch, cohort, vds_path)

    # SAMPLE QC
    sample_qc_ht_path = analysis_prefix / 'qc.ht'
    sample_qc_j = sample_qc(
        batch,
        cohort,
        vds_path,
        sample_qc_ht_path,
        depends_on=[combiner_job] if combiner_job else [],
    )

    batch.run()

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

    # def add_sample_qc_jobs(
    #     # b: hb.Batch,
    #     # mt_path: str,
    #     # samples_tsv_path: str,
    #     # sample_qc_bucket: str,
    #     # ancestry_bucket: str,
    #     # tmp_bucket: str,
    #     # analysis_bucket: str,
    #     # relatedness_bucket: str,
    #     # release_related: bool,
    #     # web_bucket: str,
    #     # filter_cutoffs_path: Optional[str],
    #     # overwrite: bool,
    #     # scatter_count: int,
    #     # combiner_job: Job,
    #     # billing_project: Optional[str] = None,
    #     # is_test: bool = False,
    #     # somalier_pairs_path: Optional[str] = None,
    #     # somalier_job: Optional[Job] = None,
    #     # highmem_workers: bool = False,
    # ):
    #     """
    #     Sub-workflow that adds sample-level QC, relatedness and ancestry analysis jobs
    #     using Hail query on the combined matrix table
    #     """
    #     job_name = 'Sample QC hard filters'
    #     hard_filtered_samples_ht_path = join(sample_qc_bucket, 'hard_filtered_samples.ht')
    #     sex_ht_path = join(sample_qc_bucket, 'sex.ht')
    #     hail_sample_qc_ht_path = join(sample_qc_bucket, 'hail_sample_qc.ht')
    #     if not can_reuse(
    #         [
    #             hard_filtered_samples_ht_path,
    #             sex_ht_path,
    #             hail_sample_qc_ht_path,
    #         ],
    #         overwrite,
    #     ):
    #         sample_qc_hardfilter_job = add_job(
    #             b,
    #             f'{utils.SCRIPTS_DIR}/sample_qc_hard_filters.py '
    #             f'--mt {mt_path} '
    #             f'--meta-tsv {samples_tsv_path} '
    #             f'{_make_filter_cutoff_param(filter_cutoffs_path, sample_qc_bucket)} '
    #             f'--out-hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
    #             f'--out-sex-ht {sex_ht_path} '
    #             f'--out-hail-sample-qc-ht {hail_sample_qc_ht_path} '
    #             f'--tmp-bucket {tmp_bucket} '
    #             + (f'--overwrite ' if overwrite else '')
    #             + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #             job_name=job_name,
    #             num_workers=scatter_count,
    #             preemptible=False,
    #             is_test=is_test,
    #             phantomjs=True,
    #             highmem=highmem_workers,
    #             depends_on=[combiner_job],
    #         )
    #     else:
    #         sample_qc_hardfilter_job = b.new_job(f'{job_name} [reuse]')
    #
    #     job_name = 'Subset MT for PCA'
    #     mt_for_pca_path = join(ancestry_bucket, 'mt_for_pca.mt')
    #     if not can_reuse(mt_for_pca_path, overwrite):
    #         subset_for_pca_job = add_job(
    #             b,
    #             f'{utils.SCRIPTS_DIR}/sample_qc_subset_mt_for_pca.py '
    #             + (f'--overwrite ' if overwrite else '')
    #             + f'--mt {mt_path} '
    #             f'--meta-tsv {samples_tsv_path} '
    #             f'--out-mt {mt_for_pca_path} '
    #             + ('--is-test ' if is_test else '')
    #             + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #             job_name=job_name,
    #             num_workers=scatter_count,
    #             is_test=is_test,
    #             phantomjs=True,
    #             highmem=highmem_workers,
    #             depends_on=[combiner_job],
    #         )
    #     else:
    #         subset_for_pca_job = b.new_job(f'{job_name} [reuse]')
    #
    #     relatedness_ht_path = join(relatedness_bucket, 'relatedness.ht')
    #     if somalier_pairs_path:
    #         job_name = 'Somalier pairs to Hail table'
    #         if not can_reuse(relatedness_ht_path, overwrite):
    #             relatedness_j = add_job(
    #                 b,
    #                 f'{utils.SCRIPTS_DIR}/sample_qc_somalier_to_ht.py '
    #                 + (f'--overwrite ' if overwrite else '')
    #                 + f'--somalier-pairs-tsv {somalier_pairs_path} '
    #                 f'--out-relatedness-ht {relatedness_ht_path} '
    #                 + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #                 job_name=job_name,
    #                 num_workers=scatter_count,
    #                 is_test=is_test,
    #                 phantomjs=True,
    #                 highmem=highmem_workers,
    #                 depends_on=[somalier_job],
    #             )
    #         else:
    #             relatedness_j = b.new_job(f'{job_name} [reuse]')
    #     else:
    #         job_name = 'Run pc_relate'
    #         if not can_reuse(relatedness_ht_path, overwrite):
    #             relatedness_ht_path = join(relatedness_bucket, 'relatedness.ht')
    #             # PC relate needs non-preemptible workes, so requesting a new cluster
    #             relatedness_j = add_job(
    #                 b,
    #                 f'{utils.SCRIPTS_DIR}/sample_qc_pcrelate.py '
    #                 + (f'--overwrite ' if overwrite else '')
    #                 + f'--pca-mt {mt_for_pca_path} '
    #                 f'--out-relatedness-ht {relatedness_ht_path} '
    #                 f'--tmp-bucket {tmp_bucket} '
    #                 + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #                 job_name=job_name,
    #                 num_workers=scatter_count,
    #                 is_test=is_test,
    #                 # Spark would have problems shuffling on preemptible workers
    #                 # (throw SparkException: Job aborted due to stage failure: ShuffleMapStage)
    #                 # so we use num_workers instead of num_secondary_workers here
    #                 preemptible=False,
    #                 highmem=highmem_workers,
    #                 depends_on=[subset_for_pca_job],
    #             )
    #         else:
    #             relatedness_j = b.new_job(f'{job_name} [reuse]')
    #
    #     cutoffs_d = utils.get_filter_cutoffs(filter_cutoffs_path)
    #
    #     job_name = 'Sample QC flag related'
    #     intermediate_related_samples_to_drop_ht_path = join(
    #         relatedness_bucket, 'intermediate_related_samples_to_drop.ht'
    #     )
    #     if not can_reuse(intermediate_related_samples_to_drop_ht_path, overwrite):
    #         flag_related_job = add_job(
    #             b,
    #             f'{utils.SCRIPTS_DIR}/sample_qc_flag_related.py '
    #             f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
    #             f'--sex-ht {sex_ht_path} '
    #             f'--relatedness-ht {relatedness_ht_path} '
    #             f'--out-ht {intermediate_related_samples_to_drop_ht_path} '
    #             f'--tmp-bucket {tmp_bucket} '
    #             f'--max-kin {cutoffs_d["pca"]["max_kin"]} '
    #             + (f'--overwrite ' if overwrite else '')
    #             + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #             job_name=job_name,
    #             num_workers=scatter_count,
    #             is_test=is_test,
    #             phantomjs=True,
    #             highmem=highmem_workers,
    #             depends_on=[sample_qc_hardfilter_job, relatedness_j],
    #         )
    #     else:
    #         flag_related_job = b.new_job(f'{job_name} [reuse]')
    #
    #     mt_for_pca_union_hgdp_path = 'gs://cpg-prophecy-main-analysis/joint-calling/v1-1/ancestry-with-hgdp/mt_for_pca_union_hgdp.mt'
    #     ancestry_bucket = join(analysis_bucket, 'ancestry_with_hgdp')
    #     num_ancestry_pcs = 16
    #     job_name = 'PCA'
    #     ancestry_analysis_bucket = ancestry_bucket
    #     ancestry_web_bucket = join(web_bucket, 'ancestry_with_hgdp')
    #     eigenvalues_ht_path = join(ancestry_analysis_bucket, f'eigenvalues.ht')
    #     scores_ht_path = join(ancestry_analysis_bucket, f'scores.ht')
    #     loadings_ht_path = join(ancestry_analysis_bucket, f'loadings.ht')
    #     inferred_pop_ht_path = join(ancestry_bucket, 'inferred_pop.ht')
    #     if not can_reuse(
    #         [eigenvalues_ht_path, scores_ht_path, loadings_ht_path, inferred_pop_ht_path],
    #         overwrite,
    #     ):
    #         pca_job = add_job(
    #             b,
    #             f'{utils.SCRIPTS_DIR}/ancestry_pca.py '
    #             f'--mt-for-pca {mt_for_pca_union_hgdp_path} '
    #             f'--meta-tsv {samples_tsv_path} '
    #             f'--min-pop-prob {cutoffs_d["pca"]["min_pop_prob"]} '
    #             f'--n-pcs {num_ancestry_pcs} '
    #             f'--out-eigenvalues-ht {eigenvalues_ht_path} '
    #             f'--out-scores-ht {scores_ht_path} '
    #             f'--out-loadings-ht {loadings_ht_path} '
    #             f'--out-inferred-pop-ht {inferred_pop_ht_path} '
    #             + (
    #                 f'--related-samples-to-drop-ht '
    #                 f'{intermediate_related_samples_to_drop_ht_path} '
    #                 if intermediate_related_samples_to_drop_ht_path
    #                 else ''
    #             )
    #             + f'--tmp-bucket {join(tmp_bucket, "ancestry")} '
    #             + (f'--overwrite ' if overwrite else '')
    #             + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #             job_name=job_name,
    #             num_workers=scatter_count,
    #             is_test=is_test,
    #             phantomjs=True,
    #             highmem=highmem_workers,
    #             depends_on=[flag_related_job, subset_for_pca_job],
    #         )
    #     else:
    #         pca_job = b.new_job(f'{job_name} [reuse]')
    #
    #     regressed_metrics_ht_path = join(sample_qc_bucket, 'regressed_metrics.ht')
    #     job_name = 'Sample QC regressed filters'
    #     if not can_reuse([regressed_metrics_ht_path], overwrite):
    #         regressed_filters_job = add_job(
    #             b,
    #             f'{utils.SCRIPTS_DIR}/sample_qc_regressed_filters.py '
    #             f'--pca-scores-ht {scores_ht_path} '
    #             f'--hail-sample-qc-ht {hail_sample_qc_ht_path} '
    #             f'--out-regressed-metrics-ht {regressed_metrics_ht_path} '
    #             f'--tmp-bucket {tmp_bucket} '
    #             + (f'--overwrite ' if overwrite else '')
    #             + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #             job_name=job_name,
    #             num_workers=scatter_count,
    #             is_test=is_test,
    #             highmem=highmem_workers,
    #             depends_on=[pca_job, sample_qc_hardfilter_job],
    #         )
    #     else:
    #         regressed_filters_job = b.new_job(f'{job_name} [reuse]')
    #
    #     job_name = 'Plot PCA and loadings'
    #     out_path_ptn = join(ancestry_web_bucket, '{scope}_pc{pci}.{ext}')
    #     paths = []
    #     for scope in ['project', 'population', 'loadings']:
    #         for ext in ['html']:
    #             paths.append(
    #                 out_path_ptn.format(scope=scope, pci=num_ancestry_pcs - 1, ext=ext)
    #             )
    #     if not can_reuse(paths, overwrite):
    #         add_job(
    #             b,
    #             f'{join(utils.SCRIPTS_DIR, "ancestry_plot.py")} '
    #             f'--meta-tsv {samples_tsv_path} '
    #             f'--eigenvalues-ht {eigenvalues_ht_path} '
    #             f'--scores-ht {scores_ht_path} '
    #             f'--loadings-ht {loadings_ht_path} '
    #             f'--inferred-pop-ht {inferred_pop_ht_path} '
    #             + f'--out-path-pattern {out_path_ptn} '
    #             + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #             job_name=job_name,
    #             num_workers=scatter_count,
    #             is_test=is_test,
    #             phantomjs=True,
    #             highmem=highmem_workers,
    #             depends_on=[pca_job, regressed_filters_job],
    #         )
    #     else:
    #         b.new_job(f'{job_name} [reuse]')
    #
    #     job_name = f'Write sample QC metadata'
    #     meta_ht_path = join(analysis_bucket, 'meta.ht')
    #     meta_tsv_path = join(analysis_bucket, 'meta.tsv')
    #     if not can_reuse(
    #         [
    #             meta_ht_path,
    #             meta_tsv_path,
    #         ],
    #         overwrite,
    #     ):
    #         metadata_qc_job = add_job(
    #             b,
    #             f'{utils.SCRIPTS_DIR}/sample_qc_write_metadata.py '
    #             f'--meta-tsv {samples_tsv_path} '
    #             f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
    #             f'--sex-ht {sex_ht_path} '
    #             f'--hail-sample-qc-ht {hail_sample_qc_ht_path} '
    #             f'--regressed-filtes-ht {regressed_metrics_ht_path} '
    #             f'--relatedness-ht {relatedness_ht_path} '
    #             + (f'--release-related ' if release_related else '')
    #             + f'--pop-ht {inferred_pop_ht_path} '
    #             f'--tmp-bucket {tmp_bucket} '
    #             f'--out-meta-ht {meta_ht_path} '
    #             f'--out-meta-tsv {meta_tsv_path} '
    #             + (f'--overwrite ' if overwrite else '')
    #             + f'{_make_filter_cutoff_param(filter_cutoffs_path, sample_qc_bucket)} '
    #             + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #             job_name=job_name,
    #             num_workers=scatter_count,
    #             is_test=is_test,
    #             highmem=highmem_workers,
    #             depends_on=[regressed_filters_job, pca_job, sample_qc_hardfilter_job],
    #         )
    #     else:
    #         metadata_qc_job = b.new_job(f'{job_name} [reuse]')
    #
    #     return metadata_qc_job, hard_filtered_samples_ht_path, meta_ht_path


if __name__ == '__main__':
    main()
