#!/usr/bin/env python

import os
from os.path import join
import setuptools


setuptools.setup(
    name='larcoh',
    version='0.1.1',
    description='Pipeline for joint calling, sample and variant QC for WGS germline '
    'variant calling data in large cohorts',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/populationgenomics/joint-calling',
    license='MIT',
    packages=['larcoh'],
    include_package_data=True,
    zip_safe=False,
    scripts=[join('scripts', fp) for fp in os.listdir('scripts') if fp.endswith('.py')],
    keywords='bioinformatics',
    install_requires=[
        'cpg-utils',
        'analysis-runner',
        'sample_metadata',
        'hail>=0.2.97',
        'cpg-gnomad',
        'pandas',
        'click',
        'google-cloud-storage',
        'google-cloud-secret-manager',
    ],
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
