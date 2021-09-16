ANALYSIS_PROJECT := tob-wgs
VERSION := v5.1
TEST_VERSION := v6-20
SCATTER_COUNT_TEST := 50
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
	--pca-pop ${PCA_POP} \
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
	--pca-pop ${PCA_POP} \
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
	--pca-pop ${PCA_POP} \
	$(REUSE_ARG)
