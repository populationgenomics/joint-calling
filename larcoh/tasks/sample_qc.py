"""
Sample-level variant QC and sex imputation
"""

from cpg_pipes.hb.batch import RegisteringBatch
from cpg_pipes.targets import Cohort
from cpg_utils import Path
from hailtop.batch.job import Job

from larcoh.pipeline_utils import dataproc_job
from larcoh.query_utils import can_reuse


def sample_qc(
    batch: RegisteringBatch,
    cohort: Cohort,
    vds_path: Path,
    out_ht_path: Path,
    depends_on: list[Job] | None = None,
) -> Job | None:
    if can_reuse(out_ht_path):
        sample_qc_hardfilter_job = None
    else:
        sample_qc_hardfilter_job = dataproc_job(
            batch,
            'sample_qc_hard_filters.py',
            params=dict(
                vds=vds_path,
                cohort_tsv=cohort.to_tsv(),
                out_ht_path=out_ht_path,
            ),
            depends_on=depends_on,
        )

    return sample_qc_hardfilter_job

    # job_name = 'Create a dense subset MT for QC'
    # subset_mt_path = analysis_prefix / 'dense-subset.mt'
    #
    # if can_reuse(subset_mt_path):
    #     subset_for_pca_job = get_batch().new_job(f'{job_name} [reuse]')
    # else:
    #     subset_for_pca_job = dataproc_job(
    #         f'sample_qc_subset_mt_for_pca.py',
    #         params=dict(
    #             mt=mt_path,
    #             cohort_tsv=cohort.to_tsv(),
    #             out_mt_path=subset_mt_path,
    #         ),
    #         depends_on=depends_on,
    #     )
    #
    # relatedness_ht_path = join(relatedness_bucket, 'relatedness.ht')
    # if somalier_pairs_path:
    #     job_name = 'Somalier pairs to Hail table'
    #     if not can_reuse(relatedness_ht_path, overwrite):
    #         relatedness_j = add_job(
    #             b,
    #             f'{utils.SCRIPTS_DIR}/sample_qc_somalier_to_ht.py '
    #             + (f'--overwrite ' if overwrite else '')
    #             + f'--somalier-pairs-tsv {somalier_pairs_path} '
    #             f'--out-relatedness-ht {relatedness_ht_path} '
    #             + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #             job_name=job_name,
    #             num_workers=scatter_count,
    #             is_test=is_test,
    #             phantomjs=True,
    #             highmem=highmem_workers,
    #             depends_on=[somalier_job],
    #         )
    #     else:
    #         relatedness_j = b.new_job(f'{job_name} [reuse]')
    # else:
    #     job_name = 'Run pc_relate'
    #     if not can_reuse(relatedness_ht_path, overwrite):
    #         relatedness_ht_path = join(relatedness_bucket, 'relatedness.ht')
    #         # PC relate needs non-preemptible workes, so requesting a new cluster
    #         relatedness_j = add_job(
    #             b,
    #             f'{utils.SCRIPTS_DIR}/sample_qc_pcrelate.py '
    #             + (f'--overwrite ' if overwrite else '')
    #             + f'--pca-mt {mt_for_pca_path} '
    #             f'--out-relatedness-ht {relatedness_ht_path} '
    #             f'--tmp-bucket {tmp_bucket} '
    #             + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #             job_name=job_name,
    #             num_workers=scatter_count,
    #             is_test=is_test,
    #             # Spark would have problems shuffling on preemptible workers
    #             # (throw SparkException: Job aborted due to stage failure: ShuffleMapStage)
    #             # so we use num_workers instead of num_secondary_workers here
    #             preemptible=False,
    #             highmem=highmem_workers,
    #             depends_on=[subset_for_pca_job],
    #         )
    #     else:
    #         relatedness_j = b.new_job(f'{job_name} [reuse]')
    #
    # cutoffs_d = utils.get_filter_cutoffs(filter_cutoffs_path)
    #
    # job_name = 'Sample QC flag related'
    # intermediate_related_samples_to_drop_ht_path = join(
    #     relatedness_bucket, 'intermediate_related_samples_to_drop.ht'
    # )
    # if not can_reuse(intermediate_related_samples_to_drop_ht_path, overwrite):
    #     flag_related_job = add_job(
    #         b,
    #         f'{utils.SCRIPTS_DIR}/sample_qc_flag_related.py '
    #         f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
    #         f'--sex-ht {sex_ht_path} '
    #         f'--relatedness-ht {relatedness_ht_path} '
    #         f'--out-ht {intermediate_related_samples_to_drop_ht_path} '
    #         f'--tmp-bucket {tmp_bucket} '
    #         f'--max-kin {cutoffs_d["pca"]["max_kin"]} '
    #         + (f'--overwrite ' if overwrite else '')
    #         + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #         job_name=job_name,
    #         num_workers=scatter_count,
    #         is_test=is_test,
    #         phantomjs=True,
    #         highmem=highmem_workers,
    #         depends_on=[sample_qc_hardfilter_job, relatedness_j],
    #     )
    # else:
    #     flag_related_job = b.new_job(f'{job_name} [reuse]')
    #
    # mt_for_pca_union_hgdp_path = 'gs://cpg-prophecy-main-analysis/joint-calling/v1-1/ancestry-with-hgdp/mt_for_pca_union_hgdp.mt'
    # ancestry_prefix = analysis_prefix / 'ancestry'
    # num_ancestry_pcs = 16
    # job_name = 'PCA'
    # ancestry_analysis_bucket = ancestry_bucket
    # ancestry_web_bucket = join(web_bucket, 'ancestry_with_hgdp')
    # eigenvalues_ht_path = join(ancestry_analysis_bucket, f'eigenvalues.ht')
    # scores_ht_path = join(ancestry_analysis_bucket, f'scores.ht')
    # loadings_ht_path = join(ancestry_analysis_bucket, f'loadings.ht')
    # inferred_pop_ht_path = join(ancestry_bucket, 'inferred_pop.ht')
    # if not can_reuse(
    #     [eigenvalues_ht_path, scores_ht_path, loadings_ht_path, inferred_pop_ht_path],
    #     overwrite,
    # ):
    #     pca_job = add_job(
    #         b,
    #         f'{utils.SCRIPTS_DIR}/ancestry_pca.py '
    #         f'--mt-for-pca {mt_for_pca_union_hgdp_path} '
    #         f'--meta-tsv {samples_tsv_path} '
    #         f'--min-pop-prob {cutoffs_d["pca"]["min_pop_prob"]} '
    #         f'--n-pcs {num_ancestry_pcs} '
    #         f'--out-eigenvalues-ht {eigenvalues_ht_path} '
    #         f'--out-scores-ht {scores_ht_path} '
    #         f'--out-loadings-ht {loadings_ht_path} '
    #         f'--out-inferred-pop-ht {inferred_pop_ht_path} '
    #         + (
    #             f'--related-samples-to-drop-ht '
    #             f'{intermediate_related_samples_to_drop_ht_path} '
    #             if intermediate_related_samples_to_drop_ht_path
    #             else ''
    #         )
    #         + f'--tmp-bucket {join(tmp_bucket, "ancestry")} '
    #         + (f'--overwrite ' if overwrite else '')
    #         + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #         job_name=job_name,
    #         num_workers=scatter_count,
    #         is_test=is_test,
    #         phantomjs=True,
    #         highmem=highmem_workers,
    #         depends_on=[flag_related_job, subset_for_pca_job],
    #     )
    # else:
    #     pca_job = b.new_job(f'{job_name} [reuse]')
    #
    # regressed_metrics_ht_path = join(sample_qc_bucket, 'regressed_metrics.ht')
    # job_name = 'Sample QC regressed filters'
    # if not can_reuse([regressed_metrics_ht_path], overwrite):
    #     regressed_filters_job = add_job(
    #         b,
    #         f'{utils.SCRIPTS_DIR}/sample_qc_regressed_filters.py '
    #         f'--pca-scores-ht {scores_ht_path} '
    #         f'--hail-sample-qc-ht {hail_sample_qc_ht_path} '
    #         f'--out-regressed-metrics-ht {regressed_metrics_ht_path} '
    #         f'--tmp-bucket {tmp_bucket} '
    #         + (f'--overwrite ' if overwrite else '')
    #         + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #         job_name=job_name,
    #         num_workers=scatter_count,
    #         is_test=is_test,
    #         highmem=highmem_workers,
    #         depends_on=[pca_job, sample_qc_hardfilter_job],
    #     )
    # else:
    #     regressed_filters_job = b.new_job(f'{job_name} [reuse]')
    #
    # job_name = 'Plot PCA and loadings'
    # out_path_ptn = join(ancestry_web_bucket, '{scope}_pc{pci}.{ext}')
    # paths = []
    # for scope in ['project', 'population', 'loadings']:
    #     for ext in ['html']:
    #         paths.append(
    #             out_path_ptn.format(scope=scope, pci=num_ancestry_pcs - 1, ext=ext)
    #         )
    # if not can_reuse(paths, overwrite):
    #     add_job(
    #         b,
    #         f'{join(utils.SCRIPTS_DIR, "ancestry_plot.py")} '
    #         f'--meta-tsv {samples_tsv_path} '
    #         f'--eigenvalues-ht {eigenvalues_ht_path} '
    #         f'--scores-ht {scores_ht_path} '
    #         f'--loadings-ht {loadings_ht_path} '
    #         f'--inferred-pop-ht {inferred_pop_ht_path} '
    #         + f'--out-path-pattern {out_path_ptn} '
    #         + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #         job_name=job_name,
    #         num_workers=scatter_count,
    #         is_test=is_test,
    #         phantomjs=True,
    #         highmem=highmem_workers,
    #         depends_on=[pca_job, regressed_filters_job],
    #     )
    # else:
    #     b.new_job(f'{job_name} [reuse]')
    #
    # job_name = f'Write sample QC metadata'
    # meta_ht_path = join(analysis_bucket, 'meta.ht')
    # meta_tsv_path = join(analysis_bucket, 'meta.tsv')
    # if not can_reuse(
    #     [
    #         meta_ht_path,
    #         meta_tsv_path,
    #     ],
    #     overwrite,
    # ):
    #     metadata_qc_job = add_job(
    #         b,
    #         f'{utils.SCRIPTS_DIR}/sample_qc_write_metadata.py '
    #         f'--meta-tsv {samples_tsv_path} '
    #         f'--hard-filtered-samples-ht {hard_filtered_samples_ht_path} '
    #         f'--sex-ht {sex_ht_path} '
    #         f'--hail-sample-qc-ht {hail_sample_qc_ht_path} '
    #         f'--regressed-filtes-ht {regressed_metrics_ht_path} '
    #         f'--relatedness-ht {relatedness_ht_path} '
    #         + (f'--release-related ' if release_related else '')
    #         + f'--pop-ht {inferred_pop_ht_path} '
    #         f'--tmp-bucket {tmp_bucket} '
    #         f'--out-meta-ht {meta_ht_path} '
    #         f'--out-meta-tsv {meta_tsv_path} '
    #         + (f'--overwrite ' if overwrite else '')
    #         + f'{_make_filter_cutoff_param(filter_cutoffs_path, sample_qc_bucket)} '
    #         + (f'--hail-billing {billing_project} ' if billing_project else ''),
    #         job_name=job_name,
    #         num_workers=scatter_count,
    #         is_test=is_test,
    #         highmem=highmem_workers,
    #         depends_on=[regressed_filters_job, pca_job, sample_qc_hardfilter_job],
    #     )
    # else:
    #     metadata_qc_job = b.new_job(f'{job_name} [reuse]')
    #
    # return sample_qc_hardfilter_jobs
