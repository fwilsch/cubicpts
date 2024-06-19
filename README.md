# Cubicpts
This repository contains a python package that provides functions to compute an exhaustive set of solutions of bounded height to several families of cubic equations as well as some utility functions to analyze them.

The functionality is mostly implemented via Cython- and C++-extensions.

## Requirements
In addition to the python dependencies, you need to have a C++-compiler supporting OpenMP and [libprimesieve](https://github.com/kimwalisch/primesieve) installed on your system and available to setuptools. In case a conda environment is active, the setup routine checks it for the presence of libprimesieve in the environment (i.e., you could satisfy the latter dependency via the conda-forge package [primesieve](https://anaconda.org/conda-forge/primesieve)).

On Macs, setuptools assumes that g++-13 is available in `/opt/homebrew` (the system-provided clang does not support OpenMP). If this is not the case on your machine, you should modify the relevant parts of `setup.py`.