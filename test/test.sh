#!/bin/bash

TIMESTAMP=$(date "+%Y%m%d-%H%M%S")

combine_gvcfs \
  --out-mt "gs://playground-au/ci/cpg_qc/${TIMESTAMP}/genomes.mt" \
  --sample-map gs://playground-au/test/gvcfs_gs.csv \
  --bucket "gs://playground-au/ci/cpg_qc/${TIMESTAMP}/" \
  --local-tmp-dir test

cpg_qc \
  --out-mt "gs://playground-au/ci/cpg_qc/${TIMESTAMP}/genomes.mt" \
  --sample-map gs://playground-au/test/gvcfs_gs.csv \
  --bucket "gs://playground-au/ci/cpg_qc/${TIMESTAMP}/" \
  --local-tmp-dir test
