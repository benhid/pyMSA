from setuptools import setup, find_packages

setup(
    name='pyMSA',
    version='1.0.0',
    description='Scoring Multiple Sequence Alignments with Python',
    url='https://github.com/benhid/pyMSA',
    author='Antonio Benítez-Hidalgo, Antonio J. Nebro',
    author_email='antonio.b@uma.es',
    license='MIT',
    python_requires='>=3',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6'
    ],
    packages=find_packages(exclude=["test.*", "tests"]),
)