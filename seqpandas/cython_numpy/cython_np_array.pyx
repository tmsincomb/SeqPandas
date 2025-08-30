# cython: language_level=3
# cython: binding=True
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

cimport cython
import numpy as np
cimport numpy as np

# Initialize NumPy C API
np.import_array()

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef to_array(list inp):
    """Convert a Python list to a NumPy int16 array using Cython."""
    cdef Py_ssize_t n = len(inp)
    cdef np.ndarray[np.int16_t, ndim=1] arr = np.zeros(n, dtype=np.int16)
    cdef Py_ssize_t idx
    
    for idx in range(n):
        arr[idx] = inp[idx]
    
    return arr
