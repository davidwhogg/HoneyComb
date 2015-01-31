import numpy
from Cython.Build import cythonize

from setuptools import setup, Extension

exts = [
    Extension("superpgram", sources=["superpgram.pyx"], include_dirs=[numpy.get_include()]), ]

setup(
    name="superpgram",
    ext_modules=cythonize(exts),
)
