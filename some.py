class Class4Tst(object):
    def __init__(self, *args, **kwargs):
        self.n = args[0]



import pyximport 
pyximport.install()
from GasLayer import foo, GasLayer
import numpy as np

if __name__ == "__main__":
    c = Class4Tst(10)
    print(np.linspace(0,1,10))
    gl = GasLayer(10)
    print(gl.get_arr())