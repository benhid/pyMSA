from os.path import abspath, dirname, join

from setuptools import find_packages, setup

basedir = abspath(dirname(__file__))

with open(join(basedir, 'README.md'), encoding='utf-8') as f:
    README = f.read()

setup(
    name='pyMSA',
    version='0.8.1',
    description='Scoring Multiple Sequence Alignments with Python',
    long_description=README,
    long_description_content_type='text/markdown',
    author='Antonio Benítez-Hidalgo, Antonio J. Nebro',
    author_email='antonio.b@uma.es',
    maintainer='Antonio Benítez-Hidalgo',
    maintainer_email='antonio.benitez@lcc.uma.es',
    license='MIT',
    url='https://github.com/benhid/pyMSA',
    packages=find_packages(exclude=["test.*", "tests"]),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Programming Language :: Python :: 3.6'
    ],
    python_requires='>=3'
)