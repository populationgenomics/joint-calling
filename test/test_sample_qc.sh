#!/bin/bash

timestamp() {
  date "+%Y%m%d-%H%M%S"
}

# Submit sample QC on the final combined matrix table
hailctl dataproc submit cpg-qc-cluster \
  --region australia-southeast1 \
  --pyfiles libs/libs.zip \
  scripts/sample_qc.py \
  --mt            gs://cpg-fewgenomes-test/cpg-qc/v2/$(timestamp)/50genomes.mt \
  --bucket        gs://cpg-fewgenomes-test/cpg-qc/v2/$(timestamp)/sample-qc \
  --out-ht        gs://cpg-fewgenomes-test/cpg-qc/v2/$(timestamp)/sample-qc.ht \
  --local-tmp-dir ~/tmp/cpg-qc/sample-qc/v2/$(timestamp) \
  --hail-billing  fewgenomes

if [ $? -eq 0 ]; then
  gsutil cp gs://cpg-fewgenomes-test/cpg-qc/v2/$(timestamp)/sample-qc.tsv .
  cut -f2 sample_qc.tsv \
    | head -n2 | tail -n1 | grep -q true || exit 1
fi
