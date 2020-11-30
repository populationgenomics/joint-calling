#!/usr/bin/env python

import os
from os.path import join
from setuptools import setup

pkg = 'qc'

try:
    import versionpy
except ImportError:
    res = input('Installation requires versionpy. Install it now? [Y/n]')
    if res.lower().startswith('n'):
        raise
    os.system('pip install versionpy')
    import versionpy

version = versionpy.get_version(pkg)
package_data = {
    pkg: versionpy.find_package_files('', pkg)
}

install_requires = []
with open("requirements.txt", "r") as requirements_file:
    for req in (line.strip() for line in requirements_file):
        if req != "hail":
            install_requires.append(req)


setup(
    name=pkg,
    script_name=pkg,
    version=version,
    description='Variant and sample QC, based on Broad\'s gnomad_qc',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url=f'https://github.com/populationgenomics/cpg_qc',
    license='MIT',
    packages=[pkg],
    package_data=package_data,
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'versionpy',
        'click',
    ],
    scripts=[
        join('scripts', 'cpg_qc.py'),
    ],
    keywords='bioinformatics',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
