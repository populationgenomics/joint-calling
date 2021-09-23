# Subset Gnomad HGDP+1KG

Takes `gs://gcp-public-data--gnomad/release/3.1/mt/genomes/gnomad.genomes.v3.1.hgdp_1kg_subset_dense.mt`, and writes into `gs://cpg-reference/hg38/ancestry`:

- `gnomad_sites_hq.mt` - subset of high quality sites for combining with a new dataset for ancestry PCA analysis
- `gnomad_sites_test_hq.mt` - same, but only 50 random samples (useful for test runs)
- `gnomad_sites_nfe_hq.mt` - only NFE samples
- `gnomad_sites_test_nfe_hq.mt` - only NFE samples out of 50 test samples

To run:

```bash
analysis-runner --dataset tob-wgs --access-level test --output-dir "ancestry/subset-gnomad/v0" --description "Subset gnomAD HGDP+1KG for PCA" python3 batch_subset_hgdp.py
```
