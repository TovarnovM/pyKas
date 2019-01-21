cdef float cfoo(float f):
    return f*f

def foo(f):
    return cfoo(f)

# import cython
import numpy as np
# cimport numpy as np 

cdef class GasLayer():
    cdef double[:] arr
    def __init__(self, n):
        self.arr = np.linspace(0,1, n)
    def get_arr(self):
        return np.asarray(self.arr)
