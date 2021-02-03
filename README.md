# Sample and variant QC

A pipeline for post-processing and filtering of population genomic variant calls.

## Installation

```sh
git clone git@github.com:populationgenomics/cpg-qc.git
pip install -e cpg-qc
```


## Usage example

```sh
cpg_qc \
  --sample-map gs://playground-au/test/gvcfs_gs.csv \
  --mt gs://playground-au/test/mt/genomes.exomesubset.mt \
  --bucket gs://playground-au/test/run \
  --local-tmp-dir test_run \
  --overwrite
```

If the matrix table specified with `--mt` doesn't exist, it will be generated automatically using [Hail's vcf_combiner](https://hail.is/docs/0.2/experimental/vcf_combiner.html).

The `--sample-map` value is a CSV file with a header as follows:

```sh
sample,population,gvcfs,contamination,alignment_summary_metrics,duplicate_metrics,insert_size_metrics,wgs_metrics
NA19238,YRI,gs://playground-au/gvcf/NA19238.g.vcf.gz,,gs://playground-au/<....>/NA19238.readgroup.alignment_summary_metrics,<...>/NA19238.duplicate_metrics,<...>/NA19238.insert_size_metrics,<...>/NA19238.wgs_metrics
```

The first column is the sample ID. The samples with data in the "population" column are used to train the random forest for population inferral of other samples.


## Description

The input is a set of GVCF files, which are merged into a Hail matrix table using vcf_combiner; then filtered on the sample level using such information as the coverage and intra-sample variant numbers/distributions; then exported into a VCF for an allele-specific VQSR, and finally exported into MT again to be filtered on the variant level using Hail.

The code is largely based on [gnomAD QC tools](https://github.com/broadinstitute/gnomad_qc), which is a collection of methods used to validate and prepare gnomAD releases. We explore gnomAD QC functions [here](docs/gnomad_qc.md). Good summaries of gnomAD QC pipeline can be found in gnomAD update blog posts:

* [https://macarthurlab.org/2017/02/27/the-genome-aggregation-database-gnomad](https://macarthurlab.org/2017/02/27/the-genome-aggregation-database-gnomad)
* [https://macarthurlab.org/2018/10/17/gnomad-v2-1](https://macarthurlab.org/2018/10/17/gnomad-v2-1)
* [https://macarthurlab.org/2019/10/16/gnomad-v3-0](https://macarthurlab.org/2019/10/16/gnomad-v3-0)
* [https://gnomad.broadinstitute.org/blog/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#sample-and-variant-quality-control](https://gnomad.broadinstitute.org/blog/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#sample-and-variant-quality-control)

Our QC pipeline:

1. gVCFs are loaded into Hail with [experimental.run_combiner](https://hail.is/docs/0.2/experimental/vcf_combiner.html)

2. Compute sample QC metrics using Hail’s [`sample_qc`](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc) module on all autosomal bi-allelic SNVs.

3. Filter outlier samples using the following cutoffs:

  * Number of SNVs: < 2.4M or > 3.75M
  * Number of singletons: > 100k
  * Hom/het ratio: > 3.3

4. Hard filtering using BAM-level metrics was performed when such metrics were available. We removed samples that were outliers for:

  * Contamination: freemix > 5% (`call-UnmappedBamToAlignedBam/UnmappedBamToAlignedBam/*/call-CheckContamination/*.selfSM`/`FREEMIX`)
  * Chimeras: > 5% (`call-AggregatedBamQC/AggregatedBamQC/*/call-CollectAggregationMetrics/*.alignment_summary_metrics`/`PCT_CHIMERAS`)
  * Duplication: > 30% (`call-UnmappedBamToAlignedBam/UnmappedBamToAlignedBam/*/call-MarkDuplicates/*.duplicate_metrics`/`PERCENT_DUPLICATION`)
  * Median insert size: < 250 (`call-AggregatedBamQC/AggregatedBamQC/*/call-CollectAggregationMetrics/*.insert_size_metrics`/`MEDIAN_INSERT_SIZE`)
  * Median coverage < 15X (`call-CollectWgsMetrics/*.wgs_metrics`/`MEDIAN_COVERAGE`)

5. Sex inferred for each sample with Hail's [`impute_sex`](https://hail.is/docs/0.2/methods/genetics.html?highlight=impute_sex#hail.methods.impute_sex). Removed samples with sex chromosome aneuploidies or ambiguous sex assignment

6. Relatedness inferred between samples using Hail's[`pc_relate`](https://hail.is/docs/0.2/methods/genetics.html?highlight=pc_relate#hail.methods.pc_relate). Identified pairs of 1st and 2nd degree relatives. Filter to a set of unrelated individuals using Hail's [`maximal_independent_set`](https://hail.is/docs/0.2/methods/misc.html?highlight=maximal_independent_set#hail.methods.maximal_independent_set) that tries to keep as many samples as possible. When multiple samples could be selected, we kept the sample with the highest coverage.

7. Going back to [`sample_qc`](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc) on all autosomal bi-allelic SNVs:

  * Filtered samples that are 4 median absolute deviations (MADs) from the median for the following metrics: `n_snp`, `r_ti_tv`, `r_insertion_deletion`, `n_insertion`, `n_deletion`, `r_het_hom_var`, `n_het`, `n_hom_var`, `n_transition`, `n_transversion`. In addition, we also filtered out samples that were in excess of 8 MADs from the median number of singletons (`n_singleton` metric).

8. For the remaining samples, export variants into a VCF.

9. Perform the allele-specific version of GATK Variant Quality Score Recalibration [VQSR](https://gatkforums.broadinstitute.org/gatk/discussion/9622/allele-specific-annotation-and-filtering), using the standard GATK training resources (HapMap, Omni, 1000 Genomes, Mills indels), with the following features:

  * SNVs:   `AS_FS`, `AS_SOR`, `AS_ReadPosRankSum`, `AS_MQRankSum`, `AS_QD`, `AS_MQ`
  * Indels: `AS_FS`, `AS_SOR`, `AS_ReadPosRankSum`, `AS_MQRankSum`, `AS_QD`

10. Apply variant hard filters:

  * No sample had a high quality genotype at this variant site (GQ>=20, DP>=10, and AB>=0.2 for heterozygotes) (all fields are populated by GATK)
  * `InbreedingCoeff` < -0.3 (there was an excess of heterozygotes at the site compared to Hardy-Weinberg expectations) (`InbreedingCoeff` is populated by GATK)
