default:
	rm -rf dist/*
	bump2version patch
	python setup.py sdist bdist_wheel
	twine upload dist/*
