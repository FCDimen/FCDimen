[![Build status](https://github.com/FCDimen/FCDimen/actions/workflows/python-app.yml/badge.svg)](https://github.com/FCDimen/FCDimen/actions)
[![PyPI](https://img.shields.io/pypi/v/fcdimen.svg)](https://pypi.python.org/pypi/fcdimen)
[![Supported Python versions](https://img.shields.io/pypi/pyversions/fcdimen.svg?logo=python&logoColor=FFE873)](https://pypi.org/project/fcdimen/)
[![Licence](https://img.shields.io/github/license/FCDimen/fcdimen.svg)](LICENSE.txt)
[![JPCL](https://img.shields.io/badge/J._Phys._Chem._Lett-Bagheri_M.%2C_Berger_E._%26_Komsa_H.P._(2023)-red)](https://doi.org/10.1021/acs.jpclett.3c01635)

# FCDimen

Python tools for analyzing dimensionality of materials structure using force constants.


## Installation

### system requirements
* Python 3.X
* Numpy >= 1.22.3
* Phonopy >= 2.8.1
* ASE >= 3.22.1
* Networkx >= 2.7.1

### Normal Installation

```bash
pip install fcdimen
```

### Developer Installation

We recommend to use a python virtual enviroment for instaling the requirements and code package:

```bash
python3 -m venv fcdimen-env
```
then activate virtual enviroment with following command:

```bash
source fcdimen-env/bin/activate
```
Get the source code from github:

```bash
git clone https://github.com/FCDimen/FCDimen.git
```

(Optional) If you want to install system requirements packages separately use following command:
```bash
cd FCDimen
pip install -r requirements.txt
```
To install the code run this command in source code directory:

```bash
pip install -e .  

```

## How to use

FCDimen needs "phonon.yaml" and "FORCE_SETS" files together or compact version of "phonon.yaml" with forces which created by phonopy.

So, First step is calculating force sets with phonopy and your favorite force calculators. You can find the full list of calculators and detailed documentation at [list of force calculators](https://phonopy.github.io/phonopy/interfaces.html).

You can find example about how to use VASP and phonopy  on [VASP & phonopy calculation](https://phonopy.github.io/phonopy/vasp.html) page.

After providing "phonon.yaml" and "FORCE_SETS" (or only compact version of "phonon.yaml") you can simply go to the directory and use fcdimen command.
You can use fcdimen -h to see more options.


## Examples

There are several examples in the examples directory that can be used like this:

```bash
fcdimen -p examples/MoS2
```
or for compact version
```bash
fcdimen -p examples/ -i MoS2.yaml
```

## How to cite

If you have used FCDimen, please cite the following article:

- "Identification of Material Dimensionality Based on Force Constant Analysis",

  Mohammad Bagheri, Ethan Berger and Hannu-Pekka Komsa, J. Phys. Chem. Lett **2023** 14 (35), 7840-7847

  https://doi.org/10.1021/acs.jpclett.3c01635  (Open Access)

  ```
  @article{FCDimen,
     author = {Bagheri, Mohammad and Berger, Ethan and Komsa, Hannu-Pekka},
     title = {Identification of Material Dimensionality Based on Force Constant Analysis}
     journal = {The Journal of Physical Chemistry Letters},
     volume = {14},
     number = {35},
     pages = {7840-7847},
     year = {2023},
     doi = {10.1021/acs.jpclett.3c01635},
     URL = {https://doi.org/10.1021/acs.jpclett.3c01635}
  }
  ```
  
## Acknowledgements

Example files are adapted from [phonondb](http://phonondb.mtl.kyoto-u.ac.jp/index.html) under CC BY 4.0.

