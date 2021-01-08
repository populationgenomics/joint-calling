#!/usr/bin/env python
from os.path import join
import setuptools

pkg = 'cpg_qc'

setuptools.setup(
    name=pkg,
    version='0.1.2',
    description='Variant and sample QC, based on Broad\'s gnomad_qc',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url=f'https://github.com/populationgenomics/{pkg}',
    license='MIT',
    packages=[pkg],
    include_package_data=True,
    zip_safe=False,
    scripts=[
        join('scripts', 'combine_gvcfs'),
        join('scripts', 'sample_qc'),
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
