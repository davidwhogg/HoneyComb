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
	       np.ndarray[DTYPE_t, ndim=1, mode="c"] data,
	       np.ndarray[DTYPE_t, ndim=1, mode="c"] ivars,
	       np.ndarray[DTYPE_t, ndim=1, mode="c"] wavenumbers):
    assert starts.shape[0] == stops.shape[0]
    assert starts.shape[0] == data.shape[0]
    assert data.shape[0] == ivars.shape[0]
    cdef int K = wavenumbers.shape[0]
    cdef int N = data.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=1, mode="c"] pgram = np.empty(K, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2, mode="c"] ATA = np.empty((3, 3), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1, mode="c"] ATy = np.empty(3, dtype=DTYPE)
    cdef int n, k
    cdef double c, s, u, d, wn, y, ivar
    for k in range(K):
        ATA[:, :] = 0.
        ATy[:] = 0.
        wn = wavenumbers[k]
        for n in range(N):
            y = data[n]
            ivar = ivars[n]
            d = starts[n]
            u = stops[n]
            c = (sin(wn * u) - sin(wn * d)) / wn
            s = (cos(wn * d) - cos(wn * u)) / wn

            ATA[0, 0] += c * c * ivar
            ATA[0, 1] += c * s * ivar
            ATA[0, 2] += c * ivar
            ATA[1, 1] += s * s * ivar
            ATA[1, 2] += s * ivar
            ATA[2, 2] += ivar
            
            ATy[0] += c * y * ivar
            ATy[1] += s * y * ivar
            ATy[2] += y * ivar

        ATA[1, 0] = ATA[0, 1]
        ATA[2, 0] = ATA[0, 2]
        ATA[2, 1] = ATA[1, 2]
    return 0.
