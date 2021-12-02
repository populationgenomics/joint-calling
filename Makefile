ANALYSIS_PROJECT := tob-wgs
VERSION := v7
TEST_VERSION := v6-26
SCATTER_COUNT_TEST := 20
SCATTER_COUNT_PROD := 50
REUSE_ARG := --reuse
# Seprately generate a specific PCA plots for this continental population
PCA_POP := nfe

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

.PHONY: 1kg_full
1kg_full:
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

.PHONY: test_to_tmp
test_to_tmp:
	analysis-runner \
	--dataset $(ANALYSIS_PROJECT) \
	--output-dir "joint-calling/test-to-tmp" \
	--description "Joint calling test-to-tmp" \
	--access-level test \
	batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_TEST) \
	--namespace tmp \
	--input-project tob-wgs \
	--analysis-project $(ANALYSIS_PROJECT) \
	--output-version ${TEST_VERSION} \
	--keep-scratch \
	--pca-pop ${PCA_POP} \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/reported_sex.tsv::1::2 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/age.csv::0::1 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::3 \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::4 \
	$(REUSE_ARG)

.PHONY: test_to_test
test_to_test:
	analysis-runner \
	--dataset $(ANALYSIS_PROJECT) \
	--output-dir "joint-calling/test-to-test" \
	--description "Joint calling test-to-test" \
	--access-level test \
	batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_TEST) \
	--namespace test \
	--input-project tob-wgs \
	--analysis-project $(ANALYSIS_PROJECT) \
	--output-version $(VERSION) \
	--keep-scratch \
	--pca-pop ${PCA_POP} \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/reported_sex.tsv::1::2 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/age.csv::0::1 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::3 \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::4 \
	$(REUSE_ARG)

.PHONY: main_to_main
main_to_main:
	analysis-runner \
	--dataset $(ANALYSIS_PROJECT) \
	--output-dir "joint-calling/main-to-main" \
	--description "Joint calling main-to-main" \
	--access-level full \
	batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_PROD) \
	--namespace main \
	--input-project tob-wgs \
	--analysis-project $(ANALYSIS_PROJECT) \
	--output-version $(VERSION) \
	--keep-scratch \
	--pca-pop ${PCA_POP} \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/reported_sex.tsv::1::2 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/age.csv::0::1 \
	--age-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::3 \
	--reported-sex-file gs://cpg-tob-wgs-main-analysis/metadata/topup_age_sex.tsv::1::4 \
	--no-filter-excesshet
	$(REUSE_ARG)

.PHONY: nagim_test
nagim_test:
	python batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_PROD) \
	--namespace test \
	--analysis-project nagim \
	--input-project thousand-genomes \
	--input-project nagim \
	--input-project amp-pd \
	--input-project hgdp \
	--input-project mgrb \
	--input-project tob-wgs \
	--input-project acute-care \
	--no-filter-excesshet \
	--source-tag nagim \
	--no-add-validation-samples \
	--keep-scratch \
	$(REUSE_ARG)
