#!/bin/bash

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

combine_gvcfs \
  --out-mt        ${DIR}/test_run/genomes.exomesubset.mt \
  --sample-map    samples.toy.csv \
  --bucket        ${DIR}/test_run/combine_gvcfs/bucket \
  --local-tmp-dir ${DIR}/test_run/combine_gvcfs/local
