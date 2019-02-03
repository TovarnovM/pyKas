# distutils: language=c++
# cython: language_level=3

cimport p1

def foo2(a):
    cdef int res = p1.cfoo(10)
    cdef p1.Bar b = p1.Bar(10,20)
    cdef double s = b.sum(7)
    return s
