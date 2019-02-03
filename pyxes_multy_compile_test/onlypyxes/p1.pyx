# distutils: language=c++
# cython: language_level=3
import numpy as np
cimport numpy as np


cdef int cfoo(int a):
    print(f"from : cfoo({a}) )")
    return a*2

cpdef int cpfoo(int a):
    print(f"from : cpfoo({a}) )")
    return cfoo(a)

def foo(a):
    print(f"from : foo({a}) )")
    return cpfoo(a)   

cdef class Bar:
    def __cinit__(self, int n, int m):
        self.a = np.zeros((n,))
        self.b = np.zeros((m,))
        cdef int i
        for i in range(n):
            self.a[i] = i*2

    cdef double sum(self, int a=10):
        cdef int i
        cdef double res = 0.0
        for i in range(a):
            res = res + self.a[i]
        return res

    def get_np(self):
        return np.array(self.a), np.array(self.b)
