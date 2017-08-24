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
included download.sh script. The package that is downloaded can include either
the entire dataset or just the raw data (i.e., not the analysis results or
figures).

```bash
# To download the entire dataset, including the results of analyses and figures
> ./download.sh
# To download only the dataset of raw data without analyses or figures
> ./download.sh raw
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
