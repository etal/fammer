#!/usr/bin/env python

"Fammer: Tools for protein superfamily sequence profile alignment and curation."

from os.path import dirname, join

DIR = (dirname(__file__) or '.')

setup_args = dict(
    name='fammer',
    version='0.2',
    description=__doc__,
    author='Eric Talevich',
    author_email='eric.talevich@gmail.com',
    url='http://github.com/etal/fammer',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=['fammerlib'],
    scripts=[join(DIR, 'fammer.py'), join(DIR, 'tmalign.py')],
)

try:
    from setuptools import setup
    setup_args.update(
        install_requires=[
            'biopython >= 1.60',
            'networkx >= 1.0',
            'biocma >= 0.1.0',
            'biofrills >= 0.2.0',
        ])
except:
    from distutils.core import setup

setup(**setup_args)

