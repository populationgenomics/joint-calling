#!/bin/bash

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

sample_qc \
  --mt            ${DIR}/test_run/genomes.exomesubset.mt \
  --meta-ht       ${DIR}/test_run/sample_qc/meta.ht \
  --sample-map    samples.toy.csv \
  --bucket        ${DIR}/test_run/sample_qc/bucket \
  --local-tmp-dir ${DIR}/test_run/sample_qc/local

cut -f2 test_run/sample_qc/meta.tsv | head -n2 | tail -n1 | grep -q true || exit 1
