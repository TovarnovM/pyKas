# distutils: language=c++
import sys
import os
cimport pyx1

def foo2(a):
    cdef int res = pyx1.foo1(33)
    return res