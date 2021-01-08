#!/bin/bash

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

sample_qc \
  --mt            ${DIR}/test_run/genomes.mt \
  --out-ht        ${DIR}/test_run/sample_qc.ht \
  --bucket        ${DIR}/test_run/sample_qc/bucket \
  --local-tmp-dir ${DIR}/test_run/sample_qc/local \
  --overwrite

if [ $? -eq 0 ]; then
  test -e ${DIR}/test_run/sample_qc.ht
  test -e ${DIR}/test_run/sample_qc.tsv
  cut -f2 ${DIR}/test_run/sample_qc.tsv \
    | head -n2 | tail -n1 | grep -q true || exit 1
fi
