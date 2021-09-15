"""
QC of newly-selected variants
"""

import os
from typing import Optional
import click
import hail as hl


GNOMAD_HGDP_1KG_MT_PATH = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)

REF_BUCKET = 'gs://cpg-reference/hg38/v1'

OUTPUT_PATH = os.path.join(
    REF_BUCKET, 'mt/gnomad.genomes.v3.1.hgdp_1kg_subset_dense{suf}.mt'
)


@click.command()
@click.option('--test', 'test', is_flag=True)
@click.option('--pop', 'pop')
def main(
    test: bool,
    pop: Optional[str],
):  # pylint: disable=missing-function-docstring
    hl.init(default_reference='GRCh38')
    mt = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT_PATH)
    suf = ''

    if test:
        suf = '_test'
        # Subset to 50 samples
        out_fpath = OUTPUT_PATH.format(suf=suf)
        if not hl.hadoop_exists(out_fpath):
            ncols = mt.count_cols()
            target_ncols = 50
            mt = mt.sample_cols(p=target_ncols / ncols, seed=42)
            mt = mt.repartition(100, shuffle=False)
            mt.write(out_fpath, overwrite=True)
        mt = hl.read_matrix_table(out_fpath)
    else:
        mt = mt.repartition(1000, shuffle=False)

    if pop:
        suf = f'{suf}_{pop}'
        out_fpath = OUTPUT_PATH.format(suf=suf)
        if not hl.hadoop_exists(out_fpath):
            # Get samples from the specified population only
            mt = mt.filter_cols(mt.population_inference.pop == pop.lower())
            mt.write(out_fpath, overwrite=True)
        mt = hl.read_matrix_table(out_fpath)

    out_fpath = OUTPUT_PATH.format(suf=f'{suf}_hq')
    if not hl.hadoop_exists(out_fpath):
        filter_high_quality_sites(
            mt,
            out_mt_path=out_fpath,
        )


def filter_high_quality_sites(
    mt: hl.MatrixTable,
    num_rows_before_ld_prune: int = 200_000,
    out_mt_path: str = None,
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

    if out_mt_path:
        mt.write(out_mt_path, overwrite=True)
        mt = hl.read_matrix_table(out_mt_path)
    return mt


if __name__ == '__main__':
    main()  # pylint: disable=E1120
