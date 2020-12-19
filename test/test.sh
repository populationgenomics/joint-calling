#!/bin/bash

TIMESTAMP=$(date "+%Y%m%d-%H%M%S")

combine_gvcfs \
  --out-mt "gs://playground-au/test/cpg_qc/combine_gvcfs/${TIMESTAMP}/genomes.exomesubset.mt" \
  --sample-map gs://playground-au/test/samples.toy.csv \
  --bucket "gs://playground-au/test/cpg_qc/combine_gvcfs/${TIMESTAMP}/work" \
  --local-tmp-dir "test_run/${TIMESTAMP}/combine_gvcfs"

sample_qc \
  --out-mt "gs://playground-au/test/cpg_qc/combine_gvcfs/${TIMESTAMP}/genomes.exomesubset.mt" \
  --sample-map gs://playground-au/test/samples.toy.csv \
  --meta-ht "gs://playground-au/test/cpg_qc/sample_qc/${TIMESTAMP}/meta.ht" \
  --bucket "gs://playground-au/test/cpg_qc/sample_qc/${TIMESTAMP}/work" \
  --local-tmp-dir "test_run/${TIMESTAMP}/sample_qc"

gsutil cp "gs://playground-au/test/cpg_qc/sample_qc/${TIMESTAMP}/meta.tsv" \
  "test_run/${TIMESTAMP}/sample_qc"

cut -f2 "test_run/${TIMESTAMP}/sample_qc/meta.tsv" | head -n2 | tail -n1 | grep -q true || exit 1

