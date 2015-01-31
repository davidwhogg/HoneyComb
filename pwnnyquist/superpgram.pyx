"""
This file is part of the HoneyComb project.
Copyright 2015 Foreman-Mackey and Hogg.

"""

from __future__ import division
cimport cython
from libc.math cimport sin, cos
import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

def superpgram(np.ndarray[DTYPE_t, ndim=1, mode="c"] starts,
               np.ndarray[DTYPE_t, ndim=1, mode="c"] stops,
	       np.ndarray[DTYPE_t, ndim=1, mode="c"] wavenumber,
	       np.ndarray[DTYPE_t, ndim=1, mode="c"] data,
	       np.ndarray[DTYPE_t, ndim=1, mode="c"] ivar):
    return 0.
