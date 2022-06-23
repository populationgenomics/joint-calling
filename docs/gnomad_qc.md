# gnomAD QC

## Loading data

### `import_vcf.py`

```sh
* import raw VCF files into Hail
      input:  ".vcf.bgz"
* left-align indel representation
* and write to MatrixTable format
      output: get_gnomad_data_path(hardcalls=False, split=True, non_refs_only=False)
          = raw_genomes_mt_path()
          = "gs://gnomad/raw/hail-*/mt/genomes/gnomad.genomes.mt"
```

### `import_resources.py`

Imports reference datasets (e.g., Clinvar annotations, methylation and CpG annotations, ExAC site annotations, truth sets) into Hail and write to MatrixTable or Table format.

```sh
input:
    "gs://gnomad-public-requester-pays/resources/grch38/purcell_5k_intervals/purcell5k.interval_list"
    "gs://gnomad-public-requester-pays/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz"
    "gs://gnomad-public-requester-pays/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed"
    "gs://gnomad-public-requester-pays/resources/grch38/syndip/full.38.20180222.vcf.gz"
    "gs://gnomad-public-requester-pays/resources/grch38/syndip/syndip.b38_20180222.bed"
    "gs://gnomad-public-requester-pays/resources/grch38/clinvar/clinvar_20190923.vcf.gz"
    "gs://gnomad-public-requester-pays/resources/grch38/dbsnp/dbsnp_b154_grch38_all_GCF_000001405.38_20200514.vcf.bgz"
    "gs://gnomad-public-requester-pays/resources/grch38/dbsnp/dbsnp_b154_grch38_all_GCF_000001405.38_20200514.vcf.header"
    "gs://gnomad-public-requester-pays/resources/grch38/dbsnp/dbsnp_b151_grch38_all_20180418.vcf.bgz"
    "gs://gnomad-public-requester-pays/resources/grch38/dbsnp/dbsnp_b151_grch38_all_20180418.vcf.header"
    "gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz"
    "gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz"
    "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    "gs://gnomad-public-requester-pays/resources/grch38/lcr_intervals/LCRFromHengHg38.txt"
    "gs://gnomad-public-requester-pays/resources/grch38/seg_dup_intervals/GRCh38_segdups.bed"
    "gs://gnomad-public-requester-pays/resources/grch38/telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.bed"

output:
    "gs://gnomad-public-requester-pays/resources/grch38/purcell_5k_intervals/purcell5k.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.mt"
    "gs://gnomad-public-requester-pays/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7_hc_regions.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/syndip/syndip.b38_20180222.mt"
    "gs://gnomad-public-requester-pays/resources/grch38/syndip/syndip_b38_20180222_hc_regions.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/clinvar/clinvar_20190923.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/dbsnp/dbsnp_b154_grch38_all_20200514.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/dbsnp/dbsnp_b151_grch38_all_20180418.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/hapmap/hapmap_3.3.hg38.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/kgp/1000G_omni2.5.hg38.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/kgp/1000G_phase1.snps.high_confidence.hg38.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/mills/Mills_and_1000G_gold_standard.indels.hg38.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/lcr_intervals/LCRFromHengHg38.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/seg_dup_intervals/GRCh38_segdups.ht"
    "gs://gnomad-public-requester-pays/resources/grch38/telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht"
```

### `load_coverage.py`

Import individual-level coverage files, write to MatrixTable and Table formats, compute summary metrics, and export summary data for release.

```sh
input:
   gnomad_meta: "gs://gnomad/metadata/genomes/gnomad.genomes.metadata.{0}.ht" -> meta_ht.bam, meta_ht.s -> bam_dict={bam->s}
      "gs://gnomad/coverage/source/genomes/coverage_files.txt":
      {base}.txt \t 'anchored' \t bam1 \t bam2 \t ...
   ["gs://gnomad/coverage/source/genomes/parts/full_{base}.gz"]
output:
   --read_coverage_files:
   ["gs://gnomad/coverage/genomes/parts/part_{base}.mt"]
   --merge_coverage_mts:
   ["gs://gnomad/coverage/hail-0.2/coverage/genomes/mt/gnomad.genomes.coverage.mt": mt.s, mt.coverage, mt.count_array, mt.locus]
   --aggregate_coverage
   "gs://gnomad-tmp/coverage/genomes_agg.ht"
   --aggregate_coverage_platforms
   "gs://gnomad/coverage/hail-0.2/coverage/genomes/mt/gnomad.genomes.coverage.groupped.mt"
   --export_coverage
   "gs://gnomad-public-requester-pays/release/2.1/coverage/genomes/gnomad.genomes.r2.1.coverage.tsv.bgz"
```

## Sample QC

### `apply_hard_filters.py`

* filter dataset to bi-allelic, high-call rate, common SNPs

```sh
input:
   get_gnomad_data(raw=True, split=False) =
      1. get_gnomad_data_path(hardcalls=False, split=False)
         = raw_genomes_mt_path()
         = "gs://gnomad/raw/hail-{0}/mt/genomes/gnomad.genomes.mt", mt.GT
      2. get_gnomad_meta(meta_version, full_meta=full_meta)
         = get_gnomad_meta_path(data_type, version)
         = metadata_genomes_ht_path()
         = "gs://gnomad/metadata/genomes/gnomad.genomes.metadata.{0}.ht"

intermediate:
   qc_mt = qc_meta_path() = '{bucket}/sample_qc/mt/gnomad.genomes.high_callrate_common_biallelic_snps{ld_pruned}.grch38.mt'
```

* import and annotate relevant metadata

```sh
input:
   qc_meta_path() -> gs://gnomad/sample_qc/input_meta/gnomad.genomes.streamlined_metadata.2018-10-10.txt.bgz, meta_ht.age
output:
   qc_mt.age
```

* infer chromosomal sex

```sh
output:
   qc_mt.is_female
   qc_mt.ambiguous_sex
```

* annotate samples failing hard filtering thresholds

```sh
input:
   qc_mt.contamination
   qc_mt.callrate
   qc_mt.pct_chimeras
   qc_mt.ambiguous_sex
   qc_mt.mean_chr20_coverage
   qc_mt.mean_dp
   qc_mt.sex_aneuploidy
   qc_mt.insert_size

   qc_mt.tcga_tumor
   qc_mt.tcga_weird_barcode
   qc_mt.tcga_below_30
   qc_mt.specific_exclusion
   qc_mt.esp
   qc_mt.non_releasable
   qc_mt.syndip

output:
   qc_mt.hard_filters
   qc_mt.perm_filters
   qc_mt.sex
   qc_mt -> qc_ht_path() = "{bucket}/sample_qc/ht/gnomad.genomeshard_filters.ht"
```

* export relevant annotations for ranking samples for removal due to relatedness (see `joint_sample_qc.py`)

```sh
input:
   qc_mt
output:
   a hardfiltered version of qc_mt, goes into rank_annotations_path() = "gs://gnomad/sample_qc/tsv/gnomad.genomes.rank_list_annotations.txt.bgz"
```

### `generate_hardcalls.py`

* generate “hardcalls” version of dataset (dropping extraneous GT fields); adjust for sex ploidy;

```sh
input:
   mt = get_gnomad_data() = raw_genomes_mt_path() = "gs://gnomad/raw/hail-{0}/mt/genomes/gnomad.genomes.mt"
   ht = qc_ht_path() = "{bucket}/sample_qc/ht/gnomad.genomeshard_filters.ht"
output:
   mt with modified GT ->
   get_gnomad_data_path(hardcalls=True, split=False) = hardcalls_mt_path(split=False) = "{bucket}/hardcalls/mt/gnomad.genomes.unsplit.mt"
```

* split multi-allelic sites;

```sh
output:
   get_gnomad_data_path(hardcalls=True, split=True)  = hardcalls_mt_path(split=False) = "{bucket}/hardcalls/mt/gnomad.genomes.mt"
```

* generate dataset version with non-reference genotypes only

```sh
output:
   get_gnomad_data_path(split=False, non_refs_only=True) = non_refs_only_mt_path(split=False) = "{bucket}/non_refs_only/mt/gnomad.genomes.unsplit.mt"
   get_gnomad_data_path(split=True, non_refs_only=True)  = non_refs_only_mt_path(split=True)  = "{bucket}/non_refs_only/mt/gnomad.genomes.mt"
```

### `joint_sample_qc.py`

* join exome and genome data;
* filter joint callset to high-call rate variants and remove rare variants;
* LD-prune joint callset;

```sh
input:
    qc_mt_path('genomes') = {bucket}/sample_qc/mt/gnomad.genomes.high_callrate_common_biallelic_snps{ld_pruned}.grch38.mt
qc_ht_path('genomes', 'hard_filters') = '{bucket}/sample_qc/ht/gnomad.genomeshard_filters.ht'
  output:
    "{bucket}/sample_qc/mt/gnomad.joint.high_callrate_common_biallelic_snps{ld_pruned}.grch38.mt"
```

* evaluate relatedness and drop related samples;

```sh
not --skip_pc_relate
output:
    "gs://gnomad/sample_qc/ht/gnomad.joint.relatedness.ht"
```

* perform genotype PCA on unrelated samples;

```sh
not --skip_relatedness
input:
    ank_annotations_path() = "gs://gnomad/sample_qc/tsv/gnomad.genomes.rank_list_annotations.txt.bgz"
    dup_pedigree_tsv_path() = "gs://gnomad/sample_qc/fam/gnomad_genomes_dup_pedigree.tsv.bgz"
output:
    rank_annotations_path() = "gs://gnomad/sample_qc/tsv/gnomad.joint.rank_list_annotations.txt.bgz"
```

* project related samples onto PCs;

```sh
input:
   compute_sample_rankings() = project_meta.ht() = "gs://gnomad/metadata/genomes_v3.1/2020-09-11_v3.1_project_meta.ht"
```

* fit random forest model on PCs for samples with known ancestry labels;
* apply RF model to assign ancestry to remaining samples;
* evaluate sample quality metrics on a joint platform-and-population basis and flag outlier samples

### `assign_subpops.py`

annotate LD-pruned joint exome and genome data with hard filters and previously inferred platform and population labels; split joint data by relatedness status; filter to samples in continental ancestry group of interest; filter relateds/unrelateds to high-call rate variants; run genotype PCA on unrelateds; project relateds onto subpop PCs; fit random forest model on subpop PCs for samples with known subpopulation labels; apply RF model to assign subpop ancestry to remaining samples

### `finalize_sample_qc.py`

combine hard filter, platform inference, relatedness, ancestry inference, subpopulation inference, and TOPMED annotations into unified metadata tables for both exomes and genomes; define release sample set; collapse small population sizes into “other” category

### `create_fam.py`

depends on relatedness output (`joint_sample_qc.py`) and hard filter flags (`apply_hard_filters.py`); infer all complete trios from kinship coefficients and sex imputation annotations, including duplicates; select representative set of unique trios in call set; evaluate inferred trios by generating random trios and comparing Mendelian violations for random trios against inferred trios; generate final list of inferred trios by filtering out trios with Mendelian violation counts exceeding standard deviation cutoff

### `get_topmed_dups.py`

depends on frequency table annotations (`annotations/generate_frequency_data.py`); create table of high-quality autosomal SNPs shared across gnomAD and TOPMED along with allele count and singleton/doubleton/tripleton annotations (on variants and samples, respectively) for identifying samples present in both gnomAD and TOPMED

## variant_qc/

### `generate_qc_annotations.py`

```sh
input:
output:
   --generate_family_stats
   {bucket}/annotations/ht/gnomad.genomes.family_stats.ht
   --generate_qc_annotations
   {bucket}/annotations/ht/gnomad.genomes.qc_stats.ht
   --generate_call_stats
   {bucket}/annotations/ht/gnomad.genomes.call_stats.ht
   --annotate_truth_data
   {bucket}/annotations/ht/gnomad.genomes.truth_data.ht
   --generate_allele_data
   {bucket}/annotations/ht/gnomad.genomes.allele_data.ht

```

### `variantqc.py`

* gather variant- and allele-level annotations used as features for training random forests model,
* impute missing entries;
* select variants for training examples;
* train random forests model using specified parameters for depth and number of trees;
* if specified, test resulting model on pre-selected region;
* save training data with metadata describing random forest parameters used;
* apply random forest model to full variant set;
* annotate results with binned scores (dependent on `create_ranked_scores.py`)
* and create variant filter status based on SNP and indel quality cutoffs, inbreeding coefficient thresholds

```sh
--annotate_for_rf
input:
   mt = get_gnomad_data() = raw_genomes_mt_path() = "gs://gnomad/raw/hail-{0}/mt/genomes/gnomad.genomes.mt"
   ht_family_stats = annotations_ht_path('family_stats') = {bucket}/annotations/ht/gnomad.genomes.family_stats.ht
   ht_qc_stats     = annotations_ht_path('qc_stats')     = {bucket}/annotations/ht/gnomad.genomes.qc_stats.ht
   ht_call_stats   = annotations_ht_path('call_stats')   = {bucket}/annotations/ht/gnomad.genomes.call_stats.ht
   call_stats_data = annotations_ht_path('call_stats')   = {bucket}/annotations/ht/gnomad.genomes.call_stats.ht
   ht_truth_data   = annotations_ht_path('truth_data')   = {bucket}/annotations/ht/gnomad.genomes.truth_data.ht
   ht_allele_data  = annotations_ht_path('allele_data')  = {bucket}/annotations/ht/gnomad.genomes.allele_data.ht
```

### `create_ranked_scores.py`

create table with rank annotations based on random forest and/or VQSR variant quality scores, with SNPs and indels handled separately, along with optional additional rankings for more specific variant types (e.g., bi-allelic variants, singletons, bi-allelic singletons); bin the variants in the rank tables by all possible rank groupings and add variant annotations from external validation datasets; compute evaluation metrics for each bin for eventual comparison of variant filtering performance across multiple random forest models and VQSR

### `calculate_concordance.py`

compute concordance against truth data (NA12878 and synthetic diploid sample) and/or for duplicate samples across exomes and genomes; create a concordance table binned by rank (both absolute and relative) for a given data type (exomes/genomes), truth sample, and filtering model (e.g., VQSR, random forest) for the purposes of evaluating variant filtering models
