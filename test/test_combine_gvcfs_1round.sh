#!/bin/bash

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

combine_gvcfs \
  --out-mt        ${DIR}/test_run/round1.mt \
  --sample-map    ${DIR}/data/sample_maps/toy-local-round1.csv \
  --bucket        ${DIR}/test_run/combine_gvcfs_round1/bucket \
  --local-tmp-dir ${DIR}/test_run/combine_gvcfs_round1/local

if [ $? -eq 0 ]; then
  test -e ${DIR}/test_run/round1.mt
  test -e ${DIR}/test_run/round1.metadata.ht
fi
