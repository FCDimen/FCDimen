# FCDimen

Python tools for analyzing dimensionality of materials structure with force constants.


## Installation

### system requirements
* Python 3.X
* Numpy
* Phonopy
* ASE
* Networkx

### Normal Installation

```bash
pip install fcdimen  
#for now:go to the package folder and  pip install .
```

### Developer Installation

We highly recommend to use a python virtual enviroment for instaling the requirements and code package. to do this you can use this command:

```bash
python3 -m venv fcdimen-env
```
then activate virtual enviroment with following command:

```bash
source fcdimen-env/bin/activate
```
Get the source code from github:

```bash
git clone https://github.com/mbagheri20/fcdimen.git
```

You can use following command for installing system requirements packages:
```bash
cd fcdimen
pip install -r requirements.txt
```
To install the code run this command in source files directory:

```bash
pip install -e .  

```

## How to use

FCDimen needs "phonon.yaml" and "FORCE_SETS" files together or compact version of "phonon.yaml" with forces which created by phonopy.
So, First step is calculating force sets with phonopy and your favorite force calculators. You can find the full list of calculators and detailed documentation at [list of force calculators](https://phonopy.github.io/phonopy/interfaces.html).
You can find example about how to use VASP and phonopy  on [VASP & phonopy calculation](https://phonopy.github.io/phonopy/vasp.html) page.
After providing "phonon.yaml" and "FORCE_SETS" (or only compact version of "phonon.yaml") you can simply go to the directory and use fcdimen command.
You can use fcdimen -h to see more options.


## examples

There are several examples in the examples directory that can be used like this:

```bash
fcdimen -p examples/MoS2
```
or for compact version
```bash
fcdimen -p examples/ -i mp-1434.yaml
```

## Acknowledgements

Example files are adapted from [phonondb](http://phonondb.mtl.kyoto-u.ac.jp/index.html) under CC BY 4.0.
