from setuptools import setup, find_packages

setup(
    name='pyMSA',
    version='0.5.1',
    description='Scoring Multiple Sequence Alignments with Python',
    author='Antonio Benítez-Hidalgo, Antonio J. Nebro',
    author_email='antonio.b@uma.es',
    maintainer='Antonio Benítez-Hidalgo',
    maintainer_email='antonio.b@uma.es',
    license='MIT',
    url='https://github.com/benhid/pyMSA',
    long_description=open('README.md').read(),
    packages=find_packages(exclude=["test.*", "tests"]),
    python_requires='>=3',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Programming Language :: Python :: 3.6'
    ]
)