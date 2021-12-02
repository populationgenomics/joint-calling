#!/usr/bin/env python

"""
Write sample QC results in 2 formats: a sample-lebel Hail Table,
and a TSV file.
"""

from os.path import join, splitext
import logging
import click
import hail as hl

from gnomad.utils.filtering import add_filters_expr

from joint_calling.utils import file_exists
from joint_calling import sample_qc as sqc, utils
from joint_calling import _version


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--meta-tsv',
    'meta_tsv_path',
    required=True,
    help='path to a TSV with QC and population metadata for the samples',
)
@click.option(
    '--hard-filtered-samples-ht',
    'hard_filtered_samples_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--hail-sample-qc-ht',
    'hail_sample_qc_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--sex-ht',
    'sex_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--custom-qc-ht',
    'custom_qc_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--regressed-filtes-ht',
    'regressed_metrics_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--relatedness-ht',
    'relatedness_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--release-related/--release-unrelated',
    'release_related',
    default=False,
    is_flag=True,
    help='Whether to keep or remove related samples from the release',
)
@click.option(
    '--pop-ht',
    'pop_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--filter-cutoffs-file',
    'filter_cutoffs_path',
    help=f'YAML file with filtering cutoffs',
)
@click.option(
    '--out-meta-ht',
    'out_meta_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=False),
)
@click.option('--out-meta-tsv', 'out_meta_tsv_path', help='Also write a TSV')
@click.option(
    '--tmp-bucket',
    'tmp_bucket',
    required=True,
    help='path to write temporary intermediate output and checkpoints',
)
@click.option(
    '--overwrite/--reuse',
    'overwrite',
    is_flag=True,
    help='if an intermediate or a final file exists, skip running the code '
    'that generates it.',
)
@click.option(
    '--hail-billing',
    'hail_billing',
    required=True,
    help='Hail billing account ID.',
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,missing-function-docstring
    meta_tsv_path: str,
    hard_filtered_samples_ht_path: str,
    hail_sample_qc_ht_path: str,
    sex_ht_path: str,
    custom_qc_ht_path: str,
    regressed_metrics_ht_path: str,
    relatedness_ht_path: str,
    release_related: bool,
    pop_ht_path: str,
    filter_cutoffs_path: str,
    out_meta_ht_path: str,
    out_meta_tsv_path: str,
    tmp_bucket: str,
    overwrite: bool,
    hail_billing: str,  # pylint: disable=unused-argument
):
    local_tmp_dir = utils.init_hail(__file__)

    cutoffs_d = utils.get_filter_cutoffs(filter_cutoffs_path)

    input_metadata_ht = utils.parse_input_metadata(meta_tsv_path, local_tmp_dir)
    hard_filtered_samples_ht = hl.read_table(hard_filtered_samples_ht_path)
    sex_ht = hl.read_table(sex_ht_path)
    custom_qc_ht = hl.read_table(custom_qc_ht_path)
    hail_sample_qc_ht = hl.read_table(hail_sample_qc_ht_path)
    regressed_metrics_ht = hl.read_table(regressed_metrics_ht_path)
    relatedness_ht = hl.read_table(relatedness_ht_path)
    pop_ht = hl.read_table(pop_ht_path)

    # Re-calculating the maximum set of unrelated samples now that
    # we have metrics adjusted for population, so newly QC-failed samples
    # are excluded
    related_samples_to_drop_ht_path = join(tmp_bucket, 'related_samples_to_drop.ht')
    related_samples_to_drop_ht = sqc.flag_related_samples(
        hard_filtered_samples_ht=hard_filtered_samples_ht,
        sex_ht=sex_ht,
        relatedness_ht=relatedness_ht,
        regressed_metrics_ht=regressed_metrics_ht,
        tmp_bucket=tmp_bucket,
        kin_threshold=cutoffs_d['pca']['max_kin'],
        out_ht_path=related_samples_to_drop_ht_path,
        overwrite=overwrite,
    )

    # Combine all intermediate tables together
    _generate_metadata(
        sample_qc_ht=hail_sample_qc_ht,
        custom_qc_ht=custom_qc_ht,
        sex_ht=sex_ht,
        input_metadata_ht=input_metadata_ht,
        hard_filtered_samples_ht=hard_filtered_samples_ht,
        regressed_metrics_ht=regressed_metrics_ht,
        pop_ht=pop_ht,
        related_samples_to_drop_after_qc_ht=related_samples_to_drop_ht,
        release_related=release_related,
        out_ht_path=out_meta_ht_path,
        out_tsv_path=out_meta_tsv_path,
        overwrite=overwrite,
    )


def _generate_metadata(
    sample_qc_ht: hl.Table,
    custom_qc_ht: hl.Table,
    sex_ht: hl.Table,
    input_metadata_ht: hl.Table,
    hard_filtered_samples_ht: hl.Table,
    regressed_metrics_ht: hl.Table,
    pop_ht: hl.Table,
    related_samples_to_drop_after_qc_ht: hl.Table,
    release_related: bool,
    out_ht_path: str,
    out_tsv_path: str,
    overwrite: bool = False,
) -> hl.Table:
    """
    Combine all intermediate tables into a single metadata table

    :param sample_qc_ht: table with a `bi_allelic_sample_qc` row field
    :param custom_qc_ht: table with a `nongnomad_snps` row field
    :param sex_ht: table with the follwing row fields:
        `f_stat`, `n_called`, `expected_homs`, `observed_homs`
    :param input_metadata_ht: table with stats from the input metadata
        (includes Picard-tools stats and assigned population tags)
    :param hard_filtered_samples_ht: table with a `hard_filters` field
    :param regressed_metrics_ht: table with global fields
        `lms`, `qc_metrics_stats`;
        and metric row fields `n_snp_residual`, `r_ti_tv_residual`, etc;
        as well as corresponding `fail_n_snp_residual`, etc,
        and a `qc_metrics_filters` row field.
    :param pop_ht: table with the following row fields:
        `pop`, `prob_CEU`, `pca_scores`, `training_pop`.
    :param related_samples_to_drop_after_qc_ht:
        table with related samples, after re-ranking them based
        on population-stratified QC
    :param release_related: whether the "release" filter should only have unrelated
        samples
    :param overwrite: overwrite checkpoints if they exist
    :return: table with relevant fields from input tables,
        annotated with the following row fields:
            'related_after_qc': bool = in `related_samples_to_drop_after_qc_ht`
            'related_before_qc': bool =
                in `related_samples_to_drop_before_qc_ht`
            'high_quality': bool =
                not `hard_filters_ht.hard_filters` &
                not `regressed_metrics_ht.qc_metrics_filters`
            'release': bool = high_quality & not `related`
    """
    logger.info('Generate the metadata HT and TSV files')

    if not overwrite and file_exists(out_ht_path):
        meta_ht = hl.read_table(out_ht_path)
    else:
        meta_ht = sex_ht.transmute(
            impute_sex_stats=hl.struct(
                f_stat=sex_ht.f_stat,
                n_called=sex_ht.n_called,
                expected_homs=sex_ht.expected_homs,
                observed_homs=sex_ht.observed_homs,
            )
        )

        meta_ht = meta_ht.annotate_globals(**regressed_metrics_ht.index_globals())

        meta_ht = meta_ht.annotate(
            sample_qc=sample_qc_ht[meta_ht.key].bi_allelic_sample_qc,
            custom_qc=custom_qc_ht[meta_ht.key],
            **hard_filtered_samples_ht[meta_ht.key],
            **regressed_metrics_ht[meta_ht.key],
            **pop_ht[meta_ht.key],
            **input_metadata_ht[meta_ht.key],
            related=hl.is_defined(related_samples_to_drop_after_qc_ht[meta_ht.key]),
        )

        meta_ht = meta_ht.annotate(
            hard_filters=hl.or_else(meta_ht.hard_filters, hl.empty_set(hl.tstr)),
        )
        if release_related:
            meta_ht = meta_ht.annotate(
                release_filters=meta_ht.hard_filters.union(meta_ht.qc_metrics_filters),
            )
        else:
            meta_ht = meta_ht.annotate(
                release_filters=add_filters_expr(
                    filters={'related': meta_ht.related},
                    current_filters=meta_ht.hard_filters.union(meta_ht.qc_metrics_filters),
                ),
            )

        meta_ht = meta_ht.annotate(
            high_quality=(hl.len(meta_ht.hard_filters) == 0),
            release=hl.len(meta_ht.release_filters) == 0,
        )
        meta_ht = meta_ht.checkpoint(out_ht_path, overwrite=True)

    out_tsv_path = out_tsv_path or splitext(out_ht_path)[0] + '.tsv'
    if overwrite or not file_exists(out_tsv_path):
        n_pcs = meta_ht.aggregate(hl.agg.min(hl.len(meta_ht.pca_scores)))
        meta_ht = meta_ht.transmute(
            **{f'PC{i + 1}': meta_ht.pca_scores[i] for i in range(n_pcs)},
            hard_filters=hl.or_missing(
                hl.len(meta_ht.hard_filters) > 0, hl.delimit(meta_ht.hard_filters)
            ),
            qc_metrics_filters=hl.or_missing(
                hl.len(meta_ht.qc_metrics_filters) > 0,
                hl.delimit(meta_ht.qc_metrics_filters),
            ),
        )
        meta_ht.flatten().export(out_tsv_path)

    return meta_ht


if __name__ == '__main__':
    main()  # pylint: disable=E1120
