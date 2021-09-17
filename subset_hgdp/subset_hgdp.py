"""
QC of newly-selected variants
"""

from os.path import join
from typing import Optional
import click
import hail as hl


GNOMAD_HGDP_1KG_MT_PATH = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)

BUCKET = 'gs://cpg-reference/hg38/ancestry/v0'
# - gnomad_sites.mt          - 105489 rows, 3942 cols
# - gnomad_sites_test.mt     - 62247 rows, 52 cols
# - gnomad_sites_test_nfe.mt - 23900 rows, 8 cols
# - gnomad_sites_nfe.mt      - 87344 rows, 675 cols
#
# Without subsetting to HQ sites:
# - all_sites/gnomad_sites_nfe.mt      - 175312130 rows, 675 cols
# - all_sites/gnomad_sites_test.mt     - 175312130 rows, 52 cols
# - all_sites/gnomad_sites_test_nfe.mt - 175312130 rows, 8 cols

NUM_TEST_SAMPLES = 50


@click.command()
@click.option('--pop', 'pop')
def main(pop: Optional[str]):  # pylint: disable=missing-function-docstring
    hl.init(default_reference='GRCh38')
    mt = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT_PATH)

    gnomad_sites_fpath = join(BUCKET, 'gnomad_sites.mt')
    if not hl.hadoop_exists(gnomad_sites_fpath):
        mt = filter_high_quality_sites(mt)
        mt = mt.repartition(5000, shuffle=False)
        mt.write(gnomad_sites_fpath, overwrite=True)
    mt = hl.read_matrix_table(gnomad_sites_fpath)

    gnomad_sites_test_fpath = join(BUCKET, 'gnomad_sites_test.mt')
    # Subset to 50 samples
    if not hl.hadoop_exists(gnomad_sites_test_fpath):
        ncols = mt.count_cols()
        target_ncols = NUM_TEST_SAMPLES
        test_mt = mt.sample_cols(p=target_ncols / ncols, seed=42)
        test_mt = test_mt.repartition(1000, shuffle=False)
        test_mt.write(gnomad_sites_test_fpath, overwrite=True)
    test_mt = hl.read_matrix_table(gnomad_sites_test_fpath)

    if pop:
        pop_fpath = join(BUCKET, f'gnomad_sites_{pop}.mt')
        if not hl.hadoop_exists(pop_fpath):
            # Get samples from the specified population only
            mt = mt.filter_cols(mt.population_inference.pop == pop.lower())
            mt = mt.repartition(1000, shuffle=False)
            mt.write(pop_fpath, overwrite=True)

        test_pop_fpath = join(BUCKET, f'gnomad_sites_test_{pop}.mt')
        if not hl.hadoop_exists(test_pop_fpath):
            test_mt = test_mt.filter_cols(
                test_mt.population_inference.pop == pop.lower()
            )
            test_mt = test_mt.repartition(1000, shuffle=False)
            test_mt.write(test_pop_fpath, overwrite=True)


def filter_high_quality_sites(
    mt: hl.MatrixTable,
    num_rows_before_ld_prune: int = 200_000,
) -> hl.MatrixTable:
    """
    1. Run `hl.variant_qc()` to calculate metrics such as AF, call rate
        and inbereeding coefficient
    2. Select variants based off of gnomAD v3 criteria: AF > 1%, call rate > 99%,
        inbreeding coefficient >-0.25 (no excess of heterozygotes)
    3. Randomly subset sites to `num_rows_before_ld_prune`
    4. LD-prune
    """
    print(f'Number of rows before filtering: {mt.count_rows()}')

    # Choose variants based off of gnomAD v3 criteria
    mt = hl.variant_qc(mt)
    mt = mt.annotate_rows(IB=hl.agg.inbreeding(mt.GT, mt.variant_qc.AF[1]))
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] > 0.01)
        & (mt.variant_qc.call_rate > 0.99)
        & (mt.IB.f_stat > -0.25)
    )

    # Randomly subsampling the matrix table to `num_rows_before_ld_prune` sites
    # before feeding it into LD prunning
    mt = mt.cache()
    nrows = mt.count_rows()
    print(f'Number of rows after filtering: {nrows}')
    if nrows > num_rows_before_ld_prune:
        print(f'Number of rows {nrows} > {num_rows_before_ld_prune}, subsetting')
        mt = mt.sample_rows(num_rows_before_ld_prune / nrows, seed=12345)

    # LD prunning
    pruned_variant_ht = hl.ld_prune(mt.GT, r2=0.1, bp_window_size=500000)
    mt = mt.filter_rows(hl.is_defined(pruned_variant_ht[mt.row_key]))
    print(f'Number of rows after prunning: {mt.count_rows()}')
    return mt


if __name__ == '__main__':
    main()  # pylint: disable=E1120
