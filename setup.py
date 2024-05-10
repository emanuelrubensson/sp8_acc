#!/usr/bin/env python
from setuptools import setup, Extension

setup(
    author="Emanuel H. Rubensson",
    name="sp8py",
    version="1.0",
    description="This is the sp8py package",
    packages=['sp8py'],
    install_requires=['numpy'],
    ext_modules=[Extension("sp8cc", sources=["sp8cc/cwrappers.cc"], libraries=['lapack'])]
)
