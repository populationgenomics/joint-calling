#!/usr/bin/env python

import os
from os.path import join, relpath
import setuptools


def find_package_files(dirpath, package, skip_exts=None):
    paths = []
    for (path, _dirs, fnames) in os.walk(join(package, dirpath)):
        for fname in fnames:
            if skip_exts and any(fname.endswith(ext) for ext in skip_exts):
                continue
            fpath = join(path, fname)
            paths.append(relpath(fpath, package))
    return paths


setuptools.setup(
    name='joint-calling',
    version='0.1.69',
    description='Pipeline for joint calling, sample and variant QC for WGS germline '
    'variant calling data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/populationgenomics/joint-calling',
    license='MIT',
    packages=['joint_calling'],
    package_data={'joint_calling': find_package_files('', 'joint_calling')},
    include_package_data=True,
    zip_safe=False,
    scripts=[
        join('scripts', 'combine_gvcfs.py'),
        join('scripts', 'sample_qc.py'),
        join('scripts', 'mt_to_vcf.py'),
        join('scripts', 'load_vqsr.py'),
        join('scripts', 'random_forest.py'),
        join('scripts', 'generate_freq_data.py'),
        join('scripts', 'generate_qc_annotations.py'),
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
