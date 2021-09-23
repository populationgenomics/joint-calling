"""
QC of newly-selected variants
"""

from os.path import join
from typing import Optional
import click
import hail as hl


LCR_INTERVALS_HT = 'gs://cpg-reference/hg38/gnomad/v0/lcr_intervals/LCRFromHengHg38.ht'
PURCELL_HT = 'gs://cpg-reference/hg38/ancestry/purcell_5k_intervals/purcell5k.ht'
GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)

BUCKET = 'gs://cpg-reference/hg38/ancestry/v1'
# - gnomad_subset.mt          - 105489 rows, 3942 cols
# - gnomad_subset_test.mt     - 62247 rows, 52 cols
# - gnomad_subset_test_nfe.mt - 23900 rows, 8 cols
# - gnomad_subset_nfe.mt      - 87344 rows, 675 cols

NUM_TEST_SAMPLES = 50


@click.command()
@click.option('--pop', 'pop')
def main(pop: Optional[str]):  # pylint: disable=missing-function-docstring
    hl.init(default_reference='GRCh38')

    gnomad_subset_mt_path = join(BUCKET, 'gnomad_subset.mt')
    if not hl.hadoop_exists(gnomad_subset_mt_path):
        ht = hl.read_table(PURCELL_HT)
        ht = ht.filter(hl.is_missing(hl.read_table(LCR_INTERVALS_HT)[ht.key]))
        mt = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
        # Filter MT to bi-allelic SNVs that are found in p5k HT
        mt = mt.filter_rows(
            (hl.len(mt.alleles) == 2)
            & hl.is_snp(mt.alleles[0], mt.alleles[1])
            & (ht[mt.locus].alleles == mt.alleles)
        )
        mt = mt.naive_coalesce(5000)
        mt.write(gnomad_subset_mt_path, overwrite=True)
    mt = hl.read_matrix_table(gnomad_subset_mt_path)

    gnomad_subset_test_fpath = join(BUCKET, 'gnomad_subset_test.mt')
    # Subset to 50 samples
    if not hl.hadoop_exists(gnomad_subset_test_fpath):
        ncols = mt.count_cols()
        target_ncols = NUM_TEST_SAMPLES
        test_mt = mt.sample_cols(p=target_ncols / ncols, seed=42)
        test_mt = test_mt.repartition(1000, shuffle=False)
        test_mt.write(gnomad_subset_test_fpath, overwrite=True)
    test_mt = hl.read_matrix_table(gnomad_subset_test_fpath)

    if pop:
        pop_fpath = join(BUCKET, f'gnomad_subset_{pop}.mt')
        if not hl.hadoop_exists(pop_fpath):
            # Get samples from the specified population only
            mt = mt.filter_cols(mt.population_inference.pop == pop.lower())
            mt = mt.repartition(1000, shuffle=False)
            mt.write(pop_fpath, overwrite=True)

        test_pop_fpath = join(BUCKET, f'gnomad_subset_test_{pop}.mt')
        if not hl.hadoop_exists(test_pop_fpath):
            test_mt = test_mt.filter_cols(
                test_mt.population_inference.pop == pop.lower()
            )
            test_mt = test_mt.repartition(1000, shuffle=False)
            test_mt.write(test_pop_fpath, overwrite=True)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
