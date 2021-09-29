"""
QC of newly-selected variants
"""

from os.path import join
from typing import Optional
import click
import hail as hl


SOMALIER_SITES_HT = 'gs://cpg-reference/hg38/somalier/v0/sites.hg38.vcf.gz'
GNOMAD_HGDP_1KG_MT = (
    'gs://gcp-public-data--gnomad/release/3.1/mt/genomes/'
    'gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt'
)


VERSION = 'v5'
BUCKET = f'gs://cpg-reference/hg38/ancestry/{VERSION}'


NUM_TEST_SAMPLES = 100

USE_SOMALIER = False


@click.command()
@click.option('--pop', 'pop')
def main(pop: Optional[str]):  # pylint: disable=missing-function-docstring
    hl.init(default_reference='GRCh38')

    sites_ht_path = join(BUCKET, 'pca_sites.ht')
    if not hl.hadoop_exists(sites_ht_path):
        if USE_SOMALIER:
            somalier_vcf_path = 'gs://cpg-reference/hg38/somalier/v0/sites.hg38.vcf.gz'
            somalier_ht = (
                hl.import_vcf(somalier_vcf_path, reference_genome='GRCh38', force=True)
                .rows()
                .key_by('locus')
            )
            somalier_ht.write(sites_ht_path)
        else:
            ht = hl.read_table('gs://cpg-reference/hg38/ancestry/v3/pca_sites_90k.ht/')
            ht.write(sites_ht_path)

    ht = hl.read_table(sites_ht_path)

    gnomad_subset_mt_path = join(BUCKET, 'gnomad_subset.mt')
    if not hl.hadoop_exists(gnomad_subset_mt_path):
        # Filter MT to bi-allelic SNVs that are found in p5k HT
        mt = hl.read_matrix_table(GNOMAD_HGDP_1KG_MT)
        mt = mt.filter_rows(
            (hl.len(mt.alleles) == 2)
            & hl.is_snp(mt.alleles[0], mt.alleles[1])
            & mt.locus.in_autosome()
            & hl.is_defined(ht[mt.locus])
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
