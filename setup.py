#coding: utf8

"""
Setup script for psgeom.
"""

from glob import glob
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(name='psgeom',
      version='0.4',
      author="TJ Lane",
      author_email="tjlane@slac.stanford.edu",
      description='scattering experiment geometry',
      packages=["psgeom"],
      package_dir={"psgeom": "psgeom"},
      scripts=[s for s in glob('scripts/*') if not s.endswith('__.py')])
