from setuptools import setup, find_packages
# import codecs
import os

# here = os.path.abspath(os.path.dirname(__file__))

# with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    # long_description = "\n" + fh.read()

VERSION = '0.0.1'
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
    long_description=long_description,
    packages=find_packages(),
    install_requires=['pymatgen', 'ase'],
    keywords=['python', 'perovskite', 'hybrid',
        'electronic structure', 'DFT', 'xTB'],
    classifiers=[
        "Development Status :: 0 - Migrating to pypi",
        "Intended Audience :: Perovskite Researchers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)
