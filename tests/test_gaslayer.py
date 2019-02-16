import sys
import os
wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
wd = os.path.dirname(os.path.dirname(wd))
sys.path.append(wd+"\\src\\")
from gaslayer import foo

def test_cimport():
    a = foo()
    assert a == 3

if __name__ == "__main__":
    for i in None:
        print(i)