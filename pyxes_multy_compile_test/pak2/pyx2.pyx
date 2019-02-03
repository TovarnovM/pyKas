import sys
import os
sys.path.append(os.path.dirname(sys.path[0])+"\\pak1")
cimport pyx1

def foo2(a):
    cdef int res = pyx1.foo1(33)
    return res