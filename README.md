# Sample and variant QC

A pipeline for post-processing and filtering of population genomic variant calls.

## Installation

```sh
git clone git@github.com:populationgenomics/joint-calling.git
cd joint-calling
conda env create -f environment-dev.yml
pip install -e .
```

## Usage with analysis-runner

You can run the workflow on test data using the [Analysis-runner](https://github.com/populationgenomics/analysis-runner).

1. Add joint-calling as submodule to your analysis repositry. E.g., to a branch "develop" of the tob-wgs project:

```sh
cd tob-wgs
git submodule add -f -b develop https://github.com/populationgenomics/joint-calling
git submodule init
git submodule update
```

1. Run the analysis runner:

```sh
# Test access level:
$ analysis-runner \
    --dataset tob-wgs \
    --output-dir "gs://cpg-tob-wgs-temporary/joint-calling-test" \
    --description "joint calling test" \
    --access-level test \
    joint-calling/workflows/drive_joint_calling.py --is-test --callset tob-wgs

# Standard access level:
$ analysis-runner \
    --dataset tob-wgs \
    --output-dir "gs://cpg-tob-wgs-temporary/joint-calling" \
    --description "joint calling prod" \
    --access-level standard \
    joint-calling/workflows/drive_joint_calling.py --callset tob-wgs \
        --version v0 --batch 0 --batch 1 --to temporary
```


## Overview

1. Find inputs. According to the specified `--dataset` and `--batch` arguments, look at `gs://cpg-<dataset>-main/gvcf/<batch-id>/` (or`gs://cpg-<dataset>-temporary/gvcf/<batch-id>/` if `--is-test`) to find GVCFs and a CSV file with QC metadata.

1. Prepare a set of GVCFs.

  - Run GATK ReblockGVCFs to annotate with allele-specific VCF INFO fields required for recalibration (QUALapprox, VarDP, RAW_MQandDP),
  - Subset GVCF to non-alt chromosomes.

1. Run the GVCF combiner using `scripts/combine_gvcfs.py`. The script interatively merges GVCFs into a sparse Matrix Table using [Hail's vcf_combiner](https://hail.is/docs/0.2/experimental/vcf_combiner.html).

1. Run `scripts/sample_qc` that performs sample-level QC, and generates a Table with filtered sample IDs, as well as a metadata Table with metrics that were used for filtering (coverage, sex, contamination, variant numbers/distributions, etc).

1. If `--run-rf` provided, use the random forest approach to perform variant QC: gather information for the random forest model, impute missing entries; select variants for training examples, train random forests model, test resulting model on pre-selected region, save training data with metadata describing random forest parameters used, and apply random forest model to the full variant set.

1. If `--run-vqsr` is specified, use the VQSR approach to perform variant QC: export the matrix table into a sites-only VCF for an allele-specific VQSR, split the VCF region-wise for parallel execution, run Gnarly Genotyper to perform "quick and dirty" joint genotyping, filter by ExcessHet, create and recalibration models for indels and SNPs separately, and evaluate the final VCF.

The sample QC and random forest variant QC pipelines are largely based on [gnomAD QC tools](https://github.com/broadinstitute/gnomad_qc), which is a collection of methods used to validate and prepare gnomAD releases. We explore gnomAD QC functions [in this document](docs/gnomad_qc.md). Good summaries of gnomAD QC pipeline can be found in gnomAD update blog posts:

- [https://macarthurlab.org/2017/02/27/the-genome-aggregation-database-gnomad](https://macarthurlab.org/2017/02/27/the-genome-aggregation-database-gnomad)
- [https://macarthurlab.org/2018/10/17/gnomad-v2-1](https://macarthurlab.org/2018/10/17/gnomad-v2-1)
- [https://macarthurlab.org/2019/10/16/gnomad-v3-0](https://macarthurlab.org/2019/10/16/gnomad-v3-0)
- [https://gnomad.broadinstitute.org/blog/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#sample-and-variant-quality-control](https://gnomad.broadinstitute.org/blog/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#sample-and-variant-quality-control)
- [https://blog.hail.is/whole-exome-and-whole-genome-sequencing-recommendations/](https://blog.hail.is/whole-exome-and-whole-genome-sequencing-recommendations/)


## Sample QC

1. GVCFs are loaded into Hail with [experimental.run_combiner](https://hail.is/docs/0.2/experimental/vcf_combiner.html)

2. Compute sample QC metrics using Hail’s [`sample_qc`](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc) module on all autosomal bi-allelic SNVs.

3. Filter outlier samples using the following cutoffs:

   - Number of SNVs: < 2.4M or > 3.75M
   - Number of singletons: > 100k
   - Hom/het ratio: > 3.3

4. Hard filtering using BAM-level metrics was performed when such metrics were available. We removed samples that were outliers for:

   - Contamination: freemix > 5% (`call-UnmappedBamToAlignedBam/UnmappedBamToAlignedBam/*/call-CheckContamination/*.selfSM`/`FREEMIX`)
   - Chimeras: > 5% (`call-AggregatedBamQC/AggregatedBamQC/*/call-CollectAggregationMetrics/*.alignment_summary_metrics`/`PCT_CHIMERAS`)
   - Duplication: > 30% (`call-UnmappedBamToAlignedBam/UnmappedBamToAlignedBam/*/call-MarkDuplicates/*.duplicate_metrics`/`PERCENT_DUPLICATION`)
   - Median insert size: < 250 (`call-AggregatedBamQC/AggregatedBamQC/*/call-CollectAggregationMetrics/*.insert_size_metrics`/`MEDIAN_INSERT_SIZE`)
   - Median coverage < 15X (`call-CollectWgsMetrics/*.wgs_metrics`/`MEDIAN_COVERAGE`)

5. Sex inferred for each sample with Hail's [`impute_sex`](https://hail.is/docs/0.2/methods/genetics.html?highlight=impute_sex#hail.methods.impute_sex). Removed samples with sex chromosome aneuploidies or ambiguous sex assignment.

6. Relatedness inferred between samples using Hail's[`pc_relate`](https://hail.is/docs/0.2/methods/genetics.html?highlight=pc_relate#hail.methods.pc_relate). Identified pairs of 1st and 2nd degree relatives. Filter to a set of unrelated individuals using Hail's [`maximal_independent_set`](https://hail.is/docs/0.2/methods/misc.html?highlight=maximal_independent_set#hail.methods.maximal_independent_set) that tries to keep as many samples as possible. When multiple samples could be selected, we kept the sample with the highest coverage.

7. PCA was a ran on high-quality variants, and RF was trained using 16 principal components as features on samples with known ancestry. Ancestry was assigned to all samples for which the probability of that ancestry was >75%.

8. [`sample_qc`](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc) was used stratified by 8 ancestry assignment PCs. Within each PC, outliers were filtered if they are 4 median absolute deviations (MADs) away from the median for the following metrics: `n_snp`, `r_ti_tv`, `r_insertion_deletion`, `n_insertion`, `n_deletion`, `r_het_hom_var`, `n_het`, `n_hom_var`, `n_transition`, `n_transversion`, or 8 MADs away from the median number of singletons (`n_singleton` metric).

9. For the remaining samples, export variants into a VCF.

10. Perform the allele-specific version of GATK Variant Quality Score Recalibration [VQSR](https://gatkforums.broadinstitute.org/gatk/discussion/9622/allele-specific-annotation-and-filtering), using the standard GATK training resources (HapMap, Omni, 1000 Genomes, Mills indels), with the following features:

   - SNVs:   `AS_FS`, `AS_SOR`, `AS_ReadPosRankSum`, `AS_MQRankSum`, `AS_QD`, `AS_MQ`
   - Indels: `AS_FS`, `AS_SOR`, `AS_ReadPosRankSum`, `AS_MQRankSum`, `AS_QD`

   - No sample had a high quality genotype at this variant site (GQ>=20, DP>=10, and AB>=0.2 for heterozygotes) (all fields are populated by GATK)
   - `InbreedingCoeff` < -0.3 (there was an excess of heterozygotes at the site compared to Hardy-Weinberg expectations) (`InbreedingCoeff` is populated by GATK)


## Development testing

For a faster iteration during development, you can also create a dataproc cluster
manually and run separate scripts on it.

```sh
gcloud config set project fewgenomes

# Start cluster
hailctl dataproc start joint-calling-cluster \
  --max-age=8h \
  --region australia-southeast1 \
  --zone australia-southeast1-a \
  --num-preemptible-workers 10 \
  --packages joint-calling,click,cpg-gnomad,google,slackclient,fsspec,sklearn,gcsfs

# Compress dependencies to upload them into the dataproc instance
# We also add gcsfs==0.3.0 to override the existing gcsfs==0.2.2
# as it required by pandas to read from the GCS. We can't specify
# gcsfs in the --packages list above as dataproc would complain
# about duplicated depenedencies. gcsfs==0.2.2 comes from
# hail/python/hailtop/hailctl/deploy.yaml
mkdir libs
cp -r joint_calling ../gnomad_methods/gnomad $CONDA_PREFIX/lib/python3.7/site-packages/gcsfs libs
cd libs
zip -r libs *
cd ..

STAMP1=$(date +"%Y-%m-%d_%H-%M-%S")

# Submit combiner job for a first set of samples
hailctl dataproc submit joint-calling-cluster \
  --region australia-southeast1 \
  --pyfiles libs/libs.zip \
  scripts/combine_gvcfs.py \
  --sample-map    gs://cpg-fewgenomes-temporary/joint-calling/50genomes-gcs-au-round1.csv \
  --out-mt        gs://cpg-fewgenomes-main/${STAMP1}/50genomes.mt \
  --bucket        gs://cpg-fewgenomes-temporary/work/vcf-combiner/${STAMP1}/ \
  --local-tmp-dir ~/tmp/joint-calling/vcf-combiner/${STAMP1}/ \
  --hail-billing  fewgenomes

STAMP2=$(date +"%Y-%m-%d_%H-%M-%S")

# Submit combiner job for a first set of samples
hailctl dataproc submit joint-calling-cluster \
  --region australia-southeast1 \
  --pyfiles libs/libs.zip \
  scripts/combine_gvcfs.py \
  --sample-map    gs://cpg-fewgenomes-temporary/joint-calling/50genomes-gcs-au-round2.csv \
  --existing-mt   gs://cpg-fewgenomes-main/${STAMP1}/50genomes.mt \
  --out-mt        gs://cpg-fewgenomes-main/${STAMP2}/50genomes.mt \
  --bucket        gs://cpg-fewgenomes-temporary/work/vcf-combiner/${STAMP2}/ \
  --local-tmp-dir ~/tmp/joint-calling/vcf-combiner/${STAMP2}/ \
  --hail-billing  fewgenomes

# Submit sample QC on the final combined matrix table
hailctl dataproc submit joint-calling-cluster \
  --region australia-southeast1 \
  --pyfiles libs/libs.zip \
  scripts/sample_qc.py \
  --mt            gs://cpg-fewgenomes-main/${STAMP2}/50genomes.mt \
  --bucket        gs://cpg-fewgenomes-main/${STAMP2}/qc/sample-qc \
  --out-ht        gs://cpg-fewgenomes-main/${STAMP2}/qc/sample-qc.ht \
  --local-tmp-dir ~/tmp/joint-calling/sample-qc/${STAMP2}/ \
  --hail-billing  fewgenomes

hailctl dataproc stop joint-calling-cluster --region australia-southeast1
```
