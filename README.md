# networkExpansionPy

Python package to construct biosphere-level metabolic networks, and run network expansion algorithms.  This package contains functions to prune biochemical reactions based on thermodynamic constraints.

The package is built from algorithms described in the following papers:

[Goldford J.E. et al, Remnants of an ancient metabolism without phosphate. Cell 168, 1–9, March 9, 2017](https://www.cell.com/fulltext/S0092-8674(17)30133-2)

[Goldford J.E. et al, Environmental boundary conditions for the origin of life converge to an organo-sulfur metabolism. Nature Ecol Evo 3,12 1715-1724, November 11,](https://pubmed.ncbi.nlm.nih.gov/31712697/)


For a guide on the use of this package see the Examples folder README.md

## Installation

In a conda or virtual environment, clone git repo and install using pip.

```sh
git clone https://github.com/jgoldford/networkExpansionPy.git
cd networkExpansionPy
pip install -e .
```
Additional steps may be necessary to fully support the package with several circumstances listed below.

For Windows with IDE:
```sh
conda install . 
```

For Windows with WSL:
```sh 
conda install conda-build
conda develop . 
```

For MacOS:
```sh
pip install . 
```
