<p align="center">
  <br/>
  <img src=resources/pymsa.png alt="pyMSA">
  <br/>
</p>

# Scoring Multiple Sequence Alignments with Python
[![Build Status](https://img.shields.io/travis/benhid/pyMSA.svg?style=flat-square)](https://travis-ci.org/benhid/pyMSA)
[![PyPI License](https://img.shields.io/pypi/l/pyMSA.svg?style=flat-square)]()
[![PyPI Python version](https://img.shields.io/pypi/pyversions/pyMSA.svg?style=flat-square)](https://pypi.org/project/pyMSA/)

pyMSA is an open source software tool aimed at providing a number of scores for
multiple sequence alignment (MSA) problems. A [tutorial](resources/tutorial-pymsa.pdf) about pyMSA is available in the resources folder of the proyect.

## Features

Score functions implemented:

* Sum of pairs,
* Star,
* Minimum entropy,
* Percentage of non-gaps,
* Percentage of totally conserved columns *and*
* STRIKE (**S**ingle s**TR**ucture **I**nduced **E**valuation).

## Downloading

To download PyMSA just clone the Git repository hosted in GitHub:
```bash
$ git clone https://github.com/benhid/pyMSA.git
$ python setup.py install
```

Alternatively, you can install it with `pip`:
```bash
$ pip install pyMSA
```

## Usage
An example of running all the included scores is located in the [`example`](examples/) folder.

<p align="center">
  <br/>
  <img src=resources/terminal.png width=600 alt="Terminal session">
  <br/>
</p>

## Authors
### Active development team
* [Antonio Benítez-Hidalgo](https://benhid.com/) <antonio.b@uma.es>
* [Antonio J. Nebro](http://www.lcc.uma.es/~antonio) <antonio@lcc.uma.es>

## License
This project is licensed under the terms of the MIT - see the [LICENSE](LICENSE) file for details.
