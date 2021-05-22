default: patch package

.PHONY: package
package:
	rm -rf dist/*
	python setup.py sdist bdist_wheel
	twine upload dist/*

.PHONY: patch
bump:
	bump2version patch

.PHONY: minor
bump:
	bump2version minor
