import numpy
from Cython.Build import cythonize

from setuptools import setup, Extension

exts = [
    Extension("superpgram._superpgram", sources=["superpgram/_superpgram.pyx"],
              include_dirs=["superpgram", numpy.get_include()]),
]

setup(
    name="superpgram",
    modules=["superpgram"],
    ext_modules=cythonize(exts),
)
