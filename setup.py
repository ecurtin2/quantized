#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=6.0', ]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

with open('requirements_dev.txt', 'r') as f:
    dev_requirements = f.readlines()

setup(
    author="Evan M Curtin",
    author_email='evanmcurtin@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Quantifying electron transit in donor-bridge-acceptor systems using probabilistic confidence.",
    entry_points={
        'console_scripts': [
            'transit_chem=transit_chem.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='transit_chem',
    name='transit_chem',
    packages=find_packages(include=['transit_chem']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    extras_require={'dev': dev_requirements},
    url='https://github.com/ecurtin2/transit_chem',
    version='0.2.0',
    zip_safe=False,
)
