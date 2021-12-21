#!/usr/bin/env python

import os
from os.path import join
import setuptools


setuptools.setup(
    name='joint-calling',
    version='0.4.20',
    description='Pipeline for joint calling, sample and variant QC for WGS germline '
    'variant calling data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/populationgenomics/joint-calling',
    license='MIT',
    packages=['joint_calling'],
    package_data={'joint_calling': ['filter_cutoffs.yaml']},
    include_package_data=True,
    zip_safe=False,
    scripts=[join('scripts', fp) for fp in os.listdir('scripts') if fp.endswith('.py')],
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
