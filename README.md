<p align="center">
  <br/>
  <img src=resources/pymsa.png alt="pyMSA">
  <br/>
</p>

# Scoring Multiple Sequence Alignments with Python
[![Build Status](https://img.shields.io/travis/benhid/pyMSA.svg?style=flat-square)](https://travis-ci.org/benhid/pyMSA)
[![PyPI License](https://img.shields.io/pypi/l/pyMSA.svg?style=flat-square)]()
[![PyPI Python version](https://img.shields.io/pypi/pyversions/pyMSA.svg?style=flat-square)]()

pyMSA is an open source software tool aimed at providing a number of scores for
multiple sequence alignment (MSA) problems.

## Features
The scores that are currently available are:
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

### STRIKE score
For computing the [STRIKE](http://www.tcoffee.org/Projects/strike/) score, 
the tool (v1.2) must be installed on the `usr/local/bin` folder.

*e.g.* After compiling run:

```bash
$ sudo cp bin/strike /usr/local/bin
```

### Substitution matrices

pyMSA has only two available substitution matrices: *PAM250*  and *Blosum62*.

Other substitution matrices can be used by reading a matrix [file](ftp://ftp.ncbi.nih.gov/blast/matrices/) with `read_matrix_from_file()` (from the `pymsa.util.substitution_matrix` module). `FileMatrix` implements this method by default:

```python
from pymsa.core.substitution_matrix import PAM250, Blosum62, FileMatrix

pam250 = PAM250()
blosum62 = Blosum62()
pam380 = FileMatrix('PAM380.txt')
```

## Usage
An example of running all the included scores is located in the [`example`](example/) folder.

## Authors
### Active development team
* Antonio Benítez-Hidalgo <antonio.b@uma.es>
* Antonio J. Nebro <antonio@lcc.uma.es>

## License
This project is licensed under the terms of the MIT - see the [LICENSE](LICENSE) file for details.
