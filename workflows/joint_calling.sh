#!/usr/bin/env bash

# To run the batch_workflow.py with the analysis runner. Assuming 
# the joint-calling repositroy is a submodule of the dataset repository.
#
# *Examples*
# 
# # Test access level - reads from `test`, writes to `temporary`:
# $ analysis-runner \
#     --dataset tob-wgs \
#     --output-dir "gs://cpg-tob-wgs-temporary/joint-calling-test" \
#     --description "joint calling test" \
#     --access-level test \
#     joint-calling/workflows/joint_calling.sh \
#     --access-level test \
#     --callset tob-wgs
# 
# Will read the inputs from:
#   gs://cpg-fewgenomes-test/gvcf/batch0/*.g.vcf.gz
# Will write the matrix tables to:
#   gs://cpg-fewgenomes-temporary/joint-vcf/v0/genomes.mt
# And other analysis files to:
#   gs://cpg-fewgenomes-temporary/joint-vcf/v0/combiner
#   gs://cpg-fewgenomes-temporary/joint-vcf/v0/random_forest
#   gs://cpg-fewgenomes-temporary/joint-vcf/v0/gvcf
# 
# # Standard access level - reads from `main`, writes to `temporary`:
# $ analysis-runner \
#     --dataset tob-wgs \
#     --output-dir "gs://cpg-tob-wgs-temporary/joint-calling" \
#     --description "joint calling standard" \
#     --access-level standard \
#     joint-calling/workflows/joint_calling.sh \
#     --access-level standard \
#     --callset tob-wgs \
#     --version v0 --batch 0 --batch 1
# 
# # Full access level - reads from `main`, writes matrix tables to `main`, 
# # other analysis files `analysis`:
# $ analysis-runner \
#     --dataset tob-wgs \
#     --output-dir "gs://cpg-tob-wgs-main/joint-calling" \
#     --description "joint calling full" \
#     --access-level full \
#     joint-calling/workflows/joint_calling.sh \
#     --access-level full \
#     --callset tob-wgs \
#     --version v0 --batch 0 --batch 1

SCRIPTPATH="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"
REPOPATH="$(dirname "$SCRIPTPATH")"

PYTHONPATH=$REPOPATH \
PATH=${REPOPATH}/scripts:$PATH \
python ${REPOPATH}/workflows/batch_workflow.py $@
