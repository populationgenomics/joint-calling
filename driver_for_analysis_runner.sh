#!/usr/bin/env bash

# To run scripts with the analysis runner. Assuming 
# the joint-calling repositroy is a submodule of the dataset repository.
# 
# Example:
# $ analysis-runner \
#     --dataset tob-wgs \
#     --output-dir "gs://cpg-tob-wgs-hail/joint-calling/test" \
#     --description "joint calling test" \
#     --access-level test \
#     joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
#     --dataset tob-wgs
#     --access-level test \

REPOPATH="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"

PYTHONPATH=$REPOPATH \
PATH=${REPOPATH}/scripts:$PATH \
python ${REPOPATH}/$@
