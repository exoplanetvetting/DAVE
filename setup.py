
#This is a minimal setup.py that I use when an OI package is missing one.
from setuptools import setup, find_packages

setup(
    name = 'Dave',
    url = 'https://github.com/exoplanetvetting/DAVE',
    version = '0.1.0',
    author = 'The Dave Team',
    author_email = 'author@gmail.com',
    description = 'Detection and Vetting of Exoplanets',

    package_dir={'': 'dave'},
)
