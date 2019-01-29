
import sys
import os
sys.path.append(os.path.dirname(sys.path[0])+"\\invariants")
cimport tube as tb

cpdef void foo():
    cdef tb.Tube t = tb.Tube([1,2,3], [2,2,3])
    print(t)
