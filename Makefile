TOB_WGS_VERSION := v7-2
TOB_WGS_TEST_VERSION := v7

REUSE_ARG := --reuse

default: patch package

.PHONY: patch
patch:
	bump2version patch
	git push

.PHONY: minor
minor:
	bump2version minor
	git push

.PHONY: package
package:
	rm -rf dist/*
	python setup.py sdist bdist_wheel
	twine upload dist/*

.PHONY: sleep
sleep:
	sleep 60

.PHONY: 1kg_main
1kg_main:
	python batch_workflow.py \
	--scatter-count 50 \
	--namespace main \
	--analysis-project thousand-genomes \
	--input-project thousand-genomes \
	--output-version v1-0 \
	--keep-scratch

.PHONY: 1kg_concordance_test
1kg_concordance_test:
	python batch_workflow.py \
	--scatter-count 20 \
	--namespace test \
	--analysis-project thousand-genomes \
	--input-project thousand-genomes \
	--output-version 1kg_concordance_test \
	--keep-scratch

.PHONY: tob_wgs_tmp
tob_wgs_tmp:
	analysis-runner \
	--dataset tob-wgs \
	--output-dir "joint-calling/test-to-tmp" \
	--description "Joint calling test-to-tmp" \
	--access-level test \
	batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_TEST) \
	--namespace tmp \
	--input-project tob-wgs \
	--analysis-project tob-wgs \
	--output-version $(TOB_WGS_TEST_VERSION) \
	--keep-scratch \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/reported_sex.tsv::1::2 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/age.csv::0::1 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::3 \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::4 \
	$(REUSE_ARG)

.PHONY: tob_wgs_test
tob_wgs_test:
	python batch_workflow.py \
	--namespace test \
	--input-project tob-wgs \
	--analysis-project tob-wgs \
	--output-version $(TOB_WGS_VERSION) \
	--keep-scratch \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/reported_sex.tsv::1::2 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/age.csv::0::1 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::3 \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::4 \
	--scatter-count 20 \
	--assume-gvcfs-are-ready \
	$(REUSE_ARG)

.PHONY: tob_wgs_main
tob_wgs_main:
	python batch_workflow.py \
	--namespace main \
	--input-project tob-wgs \
	--analysis-project tob-wgs \
	--output-version $(TOB_WGS_VERSION) \
	--keep-scratch \
	--skip-somalier \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/reported_sex.tsv::1::2 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/age.csv::0::1 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::3 \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::4 \
	--scatter-count 50 \
	--source-tag bwamem \
	--assume-gvcfs-are-ready \
	$(REUSE_ARG)

.PHONY: nagim_test
nagim_test:
	python batch_workflow.py \
	--namespace test \
	--analysis-project nagim \
	--input-project thousand-genomes \
	--input-project nagim \
	--input-project amp-pd \
	--input-project hgdp \
	--input-project mgrb \
	--input-project tob-wgs \
	--input-project acute-care \
	--source-tag nagim \
	--release-related \
	--skip-somalier \
	--no-add-validation-samples \
	--keep-scratch \
	--output-version v1-4 \
	--assume-gvcfs-are-ready \
	--scatter-count 50


# 13893 / 70 ~ 200 batches
#	--combiner-branch-factor 70 \
# processing each batch in 2 jobs
#	--combiner-batch-size 100 \ 
.PHONY: nagim_main
nagim_main:
	python batch_workflow.py \
	--namespace main \
	--analysis-project nagim \
	--input-project thousand-genomes \
	--input-project nagim \
	--input-project amp-pd \
	--input-project hgdp \
	--input-project mgrb \
	--input-project tob-wgs \
	--input-project acute-care \
	--source-tag nagim \
	--release-related \
	--skip-somalier \
	--no-add-validation-samples \
	--keep-scratch \
	--output-version v1-1 \
	--scatter-count 100 \
	--combiner-branch-factor 70 \
	--combiner-batch-size 100 \
	--highmem-workers \
	--skip-sample CPG143370 \
	--skip-sample CPG143453 \
	--skip-sample CPG143511 \
	--skip-sample CPG143743 \
	--skip-sample CPG144386 \
	--skip-sample CPG144865 \
	--skip-sample CPG145458 \
	--subset-project amp-pd \
	--subset-project mgrb
