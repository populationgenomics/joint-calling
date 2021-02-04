#!/bin/bash

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

combine_gvcfs \
  --existing-mt   ${DIR}/test_run/round1.mt \
  --out-mt        ${DIR}/test_run/genomes.mt \
  --sample-map    ${DIR}/data/sample_maps/toy-local-round2.csv \
  --bucket        ${DIR}/test_run/combine_gvcfs_round1/bucket \
  --local-tmp-dir ${DIR}/test_run/combine_gvcfs_round1/local

if [ $? -eq 0 ]; then
  test -e ${DIR}/test_run/genomes.mt
  test -e ${DIR}/test_run/genomes.metadata.ht
fi
