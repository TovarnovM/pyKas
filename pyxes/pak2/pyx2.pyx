import sys
import os
sys.path.append(os.path.dirname(sys.path[0])+"\\pak1")
from pyx1 cimport foo

def foo2(a):
    return foo(a)