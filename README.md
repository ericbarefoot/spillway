<!-- [![DOI](https://zenodo.org/badge/159710833.svg)](https://zenodo.org/badge/latestdoi/159710833) -->

# `spillway`
Discharge and water surface profiles on spillways.

`spillway.py` is module that contains the equations and solvers. `example.ipynb` is a jupyter notebook showing how to set up a model instance and use `spillway` objects.


<!-- When you use **spillway**, please cite: -->

<!-- **Wickert, A. D. and T. F. Schildgen (2019), [Long-Profile Evolution of Transport-Limited Gravel-Bed Rivers](https://www.earth-surf-dynam.net/7/17/2019/esurf-7-17-2019.html), *Earth Surf. Dynam.*, *7*, 17–43, doi:10.5194/esurf-7-17-2019.** -->

## Installation

### Via pip and PyPI

Releases will be sent to [PyPI](https://pypi.org/), and also eventually published on conda-forge.

To download and install the release version within your python system, use:

```
# Python 2
pip2 install spillway_eab # Python 2; deprecated :(

# Python 3 (implicit, assuming your packages are updated)
pip install spillway_eab

# Python 3 (explicit)
pip3 install spillway_eab
```

### Locally with pip and incorporating ongoing code modifications

To install the unreleased code from this repository and/or to make changes to it locally and have this reflected immediately in how `spillway` runs:

```
# Download the repository
git clone https://github.com/ericbarefoot/spillway.git

# Install it
# First, navigate to the root spillway directory. Then do one or both of:
pip2 install -e . # Python 2; deprecated :(
pip install -e . # Python 3; recommended
pip3 install -e . # Python 3; command to explicitly use this version
```

Of course, you may always just download the `spillway` source from here and run it as a local (rather than system-wide installed) module. But this can be inconvenient when needing to manage the directory of `spillway.py` relative to that of the driver `*.py` file that you are building to create your model run.

## Learning how to use spillway

### Flow through a spillway with constant pool levels at either end.

For a tutorial run the [Jupyter notebook contained within this package](https://github.com/ericbarefoot/spillway/example.ipynb). After installing Jupyter on your local machine, navigate to this directory in the terminal and type:
```
jupyter notebook
```
to launch it. Alternatively, a number of cloud-based services can host Jupyter notebooks.
