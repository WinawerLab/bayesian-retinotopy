# Bayesian Retinotopy
Analysis and visualization code to accompany the paper Benson and Winawer (2017).

This repository contains all the code used to analyze the data in the paper
Benson *et al*. (2017). All analysis code is written in python and can be
found in the v123 directory. The v123-analysis.ipynb iPython notebook contains
the code used to actually execute the analyses and produce the various output
files. All visualization was performed in Mathematica, and all code used to
produce the figures can be found in the visualizations.nb Mathematica
notebook. Python 2.7.9 was used to run the analyses, and Mathematica 11.0 was
used to create the visualizations.


## Authors

 * [Noah C. Benson](https://github.com/noahbenson) (*corresponding author*)
   &lt;<nben@nyu.edu>&gt;
 * [Jonathan Winawer](https://github.com/JWinawer)


## Predicting Retinotopy

If you are viewing this repository with the goal of predicting a retinotopic map
in an individual subject, you should instead look in these other places:
 * [The Neuropythy library](https://github.com/noahbenson/neuropythy); this
   Python library performs the actual Bayesian inference described in the paper.
   It can be installed in two ways; (1) you can download Neuropythy from
   GitHub, or (2) you can install Neuropythy using pip
   (`pip install neuropythy`). Neuropythy is available on
   [PyPI](https://pypi.python.org/pypi/neuropythy).
 * [The Neuropythy Docker](https://hub.docker.com/r/nben/neuropythy); this
   docker simply contains an installation of the Neuropythy library and can
   be used to easily run any of the commands included in it.

Documentation for the Neuropythy library can be found at its GitHub page; to
perform Bayesian inference with an individual subject, the required command is
the register_retinotopy command (see the Commands section of the GitHub
README file).

## Dependencies

The Python code in this repository makes heavy use of the
[Neuropythy](https://github.com/noahbenson/neuropythy) library, which should be
considered a companion of this repository.  The Mathematica notebook makes heavy
use of the [Neurotica](https://github.com/noahbenson/Neurotica) library, which
is fairly equivalent in scope to the Neuropythy library and is capable of
producing the same registration calculations. (In fact, both libraries
internally include and make use of a third [Java
library](https://github.com/noahbenson/nben) to perform the registrations.)

The versions of all libraries that were used to make the figures in the paper
are given here. The Neuropythy and Neurotica libraries used are included as
submodules of this github repository in the modules directory.
 * Neuropythy version 0.3
 * Numpy version 1.11.3
 * Scipy version 0.18.1
 * Py4j version 0.9.2

You should be able to install these dependencies automatically using setuptools;
see **Setup** below.


## Setup

### Installation

The setup.py file included with this repository should be sufficient to install
all dependencies and the 'v123' package:

```bash
> python setup.py install
```

However, this repository contains code more similar to analysis scripts than to
an actual library, so there is no particular reason to install it
permanently. An alternative is to install the Neuropythy library using pip then
to put the root of the repository on your PYTHONPATH. For example, to open the
analysis notebook:

```bash
> pip install neuropythy
> export PYTHONPATH="$PWD:$PYTHONPATH"
> jupyter-notebook-2.7 ./v123-analysis.ipynb
```

### Downloading the Data

Before running any analyses or visualization code found in either of the
included notebooks, you will have to download the data used in the paper from
the [Open Science Foundation](https://osf.io/). This can be done with the
included download.sh script. Note that the data directory is very large and
may take some time to download.

```bash
> git clone https://github.com/winawerlab/bayesian-retinotopy
# ... after cloning ...
> cd bayesian-retinotopy
# you may need to edit the permissions in order to execute download.sh
> chmod 755 ./download.sh
# To download the entire dataset
> ./download.sh
```

Note that you must have either wget or curl installed in order to use the
dowload script; if not, then a message will be printed instructing you to
download the files manually.


##License

The bayesian-retinotopy repository is free software: you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
