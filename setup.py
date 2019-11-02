#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Note: To use the 'upload' functionality of this file, you must:
#   $ pipenv install twine --dev

import io
import os
from pathlib import Path

from setuptools import find_packages, setup

# Package meta-data.
NAME = "transit_chem"
DESCRIPTION = "A Transit chem"
URL = "https://github.com/ecurtin2/transit-chem"
EMAIL = "evanmcurtin@gmail.com"
AUTHOR = "Evan Curtin"
REQUIRES_PYTHON = ">=3.4.0"
VERSION = "0.5.0"

# What packages are required for this module to be executed?
this_dir = Path(__file__).parent
REQUIRED = (this_dir / "requirements.txt").read_text().split()

# What packages are optional?
EXTRAS_DOC = (this_dir / "docs/requirements_doc.txt").read_text().split()
EXTRAS_DEV = (this_dir / "requirements_dev.txt").read_text().split()

EXTRAS = {"dev": EXTRAS_DEV, "doc": EXTRAS_DOC, "all": EXTRAS_DEV + EXTRAS_DOC}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
        long_description = "\n" + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


# Where the magic happens:
setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    # If your package is a single module, use this instead of 'packages':
    # py_modules=['mypackage'],
    # entry_points={
    #     'console_scripts': ['mycli=mymodule:cli'],
    # },
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license="MIT",
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],
)
