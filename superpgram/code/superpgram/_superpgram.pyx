"""
This file is part of the HoneyComb project.
Copyright 2015 Foreman-Mackey and Hogg.
"""

from __future__ import division
cimport cython
from libc.math cimport sin, cos
import numpy as np
cimport numpy as np
import scipy.linalg.lapack

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

# Pointer to LAPACK function.
cdef extern from "f2pyptr.h":
    void *f2py_pointer(object) except NULL
ctypedef int dposv_t(
        char* uplo, int* n, int* nrhs,
        double* a, int* lda,
        double* b, int* ldb,
        int* info)
cdef dposv_t* dposv = <dposv_t*>f2py_pointer(scipy.linalg.lapack.dposv._cpointer)


@cython.boundscheck(False)
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
    cdef np.ndarray[DTYPE_t, ndim=2, mode="fortran"] ATA = \
            np.empty((3, 3), dtype=DTYPE, order="F")
    cdef np.ndarray[DTYPE_t, ndim=1, mode="fortran"] ATy = \
            np.empty(3, dtype=DTYPE, order="F")
    cdef int n, k
    cdef double c, s, u, d, wn, y, ivar

    cdef int info = 0
    cdef int one = 1
    cdef int three = 3

    for k in range(K):
        ATA[:, :] = 0.
        ATy[:] = 0.
        wn = wavenumbers[k]
        for n in range(N):
            y = data[n]
            ivar = ivars[n]
            d = starts[n]
            u = stops[n]
            # go trig identity
            meanangle = 0.5 * wn * (u + d)
            sortadeltat = 2. * sin(0.5 * wn * (u - d)) / wn
            c = cos(meanangle) * sortadeltat
            s = sin(meanangle) * sortadeltat
            # the above four lines replace the following two lines
            # c = (sin(wn * u) - sin(wn * d)) / wn
            # s = (cos(wn * d) - cos(wn * u)) / wn

            ATA[0, 0] += c * c * ivar
            ATA[0, 1] += c * s * ivar
            ATA[0, 2] += c * ivar
            ATA[1, 1] += s * s * ivar
            ATA[1, 2] += s * ivar
            ATA[2, 2] += ivar

            ATy[0] += c * y * ivar
            ATy[1] += s * y * ivar
            ATy[2] += y * ivar

        dposv("U", &three, &one, <double*>ATA.data, &three,
              <double*>ATy.data, &three, &info)

        if info:
            raise np.linalg.LinAlgError("Solve failed with {0}".format(info))

        pgram[k] = ATy[0] * ATy[0] + ATy[1] * ATy[1]
    return pgram
