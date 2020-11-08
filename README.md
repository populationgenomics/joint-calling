# CPG QC

Pipeline for post-processing the population genomics samples and variants. The input is a set of gVCFs, which are merged into a Hail matrix table using vcf_combiner; then filtered on the sample level using such information as the coverage and intra-sample variant numbers/distributions; then exported into a VCF for an allele-specific VQSR, and finally exported into MT again to be filtered on the variant level using Hail.

The code is largely based on [gnomAD QC tools](https://github.com/broadinstitute/gnomad_qc), which is a collection of methods used to validate and prepare gnomAD releases. We [summarize](gnomad_qc.md) gnomAD QC here.


