<p align="center">
  <br/>
  <img src=resources/pymsa.png alt="PyMSA">
  <br/>
</p>

# Scoring Multiple Sequence Alignments with Python
[![Build Status](https://travis-ci.org/benhid/pyMSA.svg?branch=master)](https://travis-ci.org/benhid/pyMSA)
[![PyPI](https://img.shields.io/pypi/l/pyMSA.svg)]()
[![PyPI](https://img.shields.io/pypi/v/pyMSA.svg)]()

PyMSA is an open source software tool aimed at providing a number of scores for
multiple sequence alignment (MSA) problems.

## Features
The scores that are currently available are:
* Sum of pairs,
* Star,
* Minimum entropy,
* Percentage of non-gaps,
* Percentage of totally conserved columns *and*
* [STRIKE](http://www.tcoffee.org/Projects/strike/) (*S*ingle s*TR*ucture *I*nduced *E*valuation).

## Downloading
To download PyMSA just clone the Git repository hosted in GitHub:
```bash
$ git clone https://github.com/benhid/pyMSA.git
$ python setup.py install
```

Alternatively, you can install with pip:
```bash
$ pip install pyMSA
```

## Requirements
PyMSA has been developed with Python 3.6.0 :: [Anaconda](https://www.continuum.io) 4.3.0 (x86_64).

To install all dependencies use:
```bash
$ pip install -r requirements.txt
```

## Usage
An example of running all the included scores is located in the [`example`](example/) folder.

## Authors
### Active development team
* Antonio Benítez <antonio.b@uma.es>
* Antonio J. Nebro <antonio@lcc.uma.es>

## License
This project is licensed under the terms of the MIT - see the [LICENSE](LICENSE) file for details.
