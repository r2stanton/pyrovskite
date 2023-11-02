from setuptools import setup, find_packages
import os

VERSION = '1.0.0'
DESCRIPTION = 'Python package for 2D- and 3D-perovskites'
LONG_DESCRIPTION = 'A software package for the high throughput construction, analysis, and featurization of two- and three-dimensional perovskite systems.'

# Setting up
setup(
    name="pyrovskite",
    version=VERSION,
    author="Robert Stanton, Dhara J. Trivedi",
    author_email="<stantor@clarkson.edu>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=['pymatgen', 'ase'],
    keywords=['python', 'perovskite', 'hybrid',
        'electronic structure', 'DFT', 'xTB'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)
