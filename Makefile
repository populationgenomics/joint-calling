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
