# scgenome: single cell whole genome analysis in python

scgenome is scalable python toolkit for analyzing single-cell whole genome
data built on [anndata](https://anndata.readthedocs.io) and inspired by
[scanpy](https://scanpy.readthedocs.io).  scgenome includes preprocessing,
visualization, and clustering functionality and can be used to analyze
copy number, allele specific copy number, SNV and breakpoint features.

## Installation

Installation from pypi is recommended:

```
pip install scgenome
```

We also recommend installing scgenome into a virtual environment using
virtualenv.  To create a fresh environment with scgenome use the following
steps:

```
virtualenv venv
source venv/bin/activate
pip install scgenome
```

To install from source, clone this repo and use the following steps in
the repo directory:

```
virtualenv venv
source venv/bin/activate
python setup.py develop
```

## Documentation

Coming Soon

## Pip build

To build with pip and distribute to pypi, use the following commands:

    python setup.py build_ext --force sdist
    twine upload --repository pypi dist/*


