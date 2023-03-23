from setuptools import setup
import py2exe

setup(
    name='MEM4_outliers',
    version='',
    packages=['pqtrees', 'pqtrees.tests', 'pqtrees.tests.pqtree', 'pqtrees.utilities', 'pqtrees.pqtree_helpers',
              'pqtrees.common_intervals'],
    url='',
    license='',
    author='Eden',
    author_email='edenozery@gmail.com',
    description='',
    console=['main.py']
)
