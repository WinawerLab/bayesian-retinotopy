#! /usr/bin/env python

import os
from setuptools import setup

version = '0.0.1'

setup(
    name='v123-retinotopic-maps',
    version=version,
    description='Analysis repository that is the companion of the paper Benson et al. (2017)',
    keywords='neuroscience mesh cortex registration',
    author='Noah C. Benson',
    author_email='nben@nyu.edu',
    maintainer_email='nben@nyu.edu',
    long_description='''
                     See the README.md file at the github repository for this package:
                     https://github.com/winawerlab/v123-retinotopic-maps
                     ''',
    url='https://github.com/winawerlab/v123-retinotopic-maps',
    download_url='https://github.com/winawerlab/v123-retinotopic-maps',
    license='GPLv3',
    classifiers=['Intended Audience :: Science/Research',
                 'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                 'Programming Language :: Python :: 2.7',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Information Analysis',
                 'Topic :: Scientific/Engineering :: Medical Science Apps.',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Operating System :: Unix',
                 'Operating System :: MacOS'],
    packages=['v123'],
    include_package_data=True,
    package_data={'': ['LICENSE.txt']},
    install_requires=['neuropythy>=0.2',
                      'numpy>=1.2',
                      'scipy>=0.7'])
