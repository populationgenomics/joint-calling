VERSION := v5.1
TEST_VERSION := v6-14
SCATTER_COUNT_TEST := 50
SCATTER_COUNT_PROD := 50
ANALYSIS_PROJECT := tob-wgs
REUSE_ARG := --reuse
PRE_COMPUTED_HGDP := gs://cpg-tob-wgs-test-tmp/joint-calling/v6-11/sample_qc/mt_subset_for_pca_with_hgdp.mt

default: patch package

.PHONY: patch
patch:
	bump2version patch

.PHONY: minor
minor:
	bump2version minor

.PHONY: package
package:
	rm -rf dist/*
	python setup.py sdist bdist_wheel
	twine upload dist/*

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
	--analysis-project $(ANALYSIS_PROJECT) \
	--input-project tob-wgs \
	--output-version ${TEST_VERSION} \
	--keep-scratch \
	--pre-computed-hgdp-union-mt $(PRE_COMPUTED_HGDP) \
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
	--input-project tob-wgs \
	--output-version $(VERSION) \
	--keep-scratch \
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
	--input-project tob-wgs \
	--output-version $(VERSION) \
	$(REUSE_ARG)
