from setuptools import setup, find_packages

setup(
    name='pyMSA',
    version='0.1.1',
    description='Scoring Multiple Sequence Alignments with Python',
    url='https://github.com/benhid/pyMSA',
    author='Antonio Benitez, Antonio J. Nebro',
    author_email='antonio.b@uma.es',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6'
    ],
    packages=find_packages(exclude=['*Test']),
)