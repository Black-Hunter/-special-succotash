# compiler directives for cython compiler
# distutils: language = c++
# distutils: sources = ../svd.cpp
# cython: language_level = 3
# cython: boundscheck = False
# cython: wraparound = False

# import numpy (we want to return numpy arrays)
import numpy as np

# we tell cython what from the header "../svd.h" we want to use
cdef extern from "svd.h":
    struct SVD:
        const double* U  # m * k
        const double* s  # k
        const double* VT  # k * n
        const int k
        const int m
        const int n

    void free_svd_result(SVD * r)

    SVD * compute_svd(double* A, int m, int n,
                            double threshold, double accuracy)

def partial_svd(double[:,::1] A, double threshold=-1.0, double accuracy=0.00000001):
    # solve svd
    cdef SVD * result = compute_svd(&A[0,0], A.shape[0],
                                           A.shape[1], threshold, accuracy)
    # number of rows, columns in result matrices
    cdef int m = result[0].m, k = result[0].k, n = result[0].n
	
    # create empty numpy arrays
    VT = np.empty(shape=(n, k), dtype=np.double)
    s  = np.empty(shape=(k), dtype=np.double)
    U  = np.empty(shape=(k, m), dtype=np.double)

    # copy results into numpy arrays
    # could be done without copying, but code would be complicated!
    if(k > 0):
        U[:,:] = <double[:k,:m]> result[0].U
        s[:] = <double[:k]> result[0].s
        VT[:,:] = <double[:n,:k]> result[0].VT
    
    # free memory
    free_svd_result(result)
    
    # return the three numpy arrays
    return U, s, VT
