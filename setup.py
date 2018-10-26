# setup.py

#from distutils.core import setup
from setuptools import setup, find_packages

setup(name = 'mylib',
      version = '1.0',
      packages = ['mylib'],
      package_dir = {'mylib':'mylib'},
      include_package_data = True
     )
