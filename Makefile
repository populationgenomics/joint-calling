VERSION := v5
TEST_VERSION := v6-4
SCATTER_COUNT_TEST := 10
SCATTER_COUNT_PROD := 100
ANALYSIS_PROJECT := tob-wgs
REUSE_ARG := --reuse

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
	--description "Joint calling test-to-temporary" \
	--access-level test \
	workflows/batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_TEST) \
	--namespace tmp \
	--batch batch1 \
	--analysis-project $(ANALYSIS_PROJECT) \
	--input-project tob-wgs-test \
	--output-version ${TEST_VERSION} \
	--keep-scratch \
	$(REUSE_ARG)

.PHONY: test_to_test
test_to_test:
	analysis-runner \
	--dataset $(ANALYSIS_PROJECT) \
	--output-dir "joint-calling/test-to-test" \
	--description "Joint calling test-to-test" \
	--access-level test \
	workflows/batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_TEST) \
	--namespace test \
	--input-project tob-wgs-test \
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
	workflows/batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_PROD) \
	--namespace main \
	--input-project tob-wgs \
	--analysis-project $(ANALYSIS_PROJECT) \
	--input-project tob-wgs \
	--output-version $(VERSION) \
	$(REUSE_ARG)
