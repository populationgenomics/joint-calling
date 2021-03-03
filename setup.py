#!/usr/bin/env python

from os.path import join
import setuptools

PKG = 'cpg_qc'

setuptools.setup(
    name=PKG,
    version='0.1.2',
    description='Variant and sample QC, based on Broad\'s gnomad_qc',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url=f'https://github.com/populationgenomics/{PKG}',
    license='MIT',
    packages=[PKG],
    include_package_data=True,
    zip_safe=False,
    scripts=[
        join('scripts', 'combine_gvcfs.py'),
        join('scripts', 'sample_qc.py'),
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
