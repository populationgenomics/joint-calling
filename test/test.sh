#!/bin/bash

cpg_qc \
  --mt gs://playground-au/test/mt/genomes.exomesubset.mt \
  --sample-map gs://playground-au/test/gvcfs_gs.csv \
  --bucket gs://playground-au/test/run \
  --local-tmp-dir test \
  --overwrite
