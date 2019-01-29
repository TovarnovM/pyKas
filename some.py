from invariants.tube import Tube
import numpy as np
import sys

t = Tube([1,2,3], [1,1,1])
print(np.asarray(t.get_xs()))
print(sys.path[0])

# class Class4Tst(object):
#     def __init__(self, *args, **kwargs):
#         self.n = args[0]

#     def foo(self, b):
#         return b

    



# import pyximport 
# pyximport.install()
# from GasLayer import foo, GasLayer
# import numpy as np

# if __name__ == "__main__":
#     c = Class4Tst(10)
#     print(np.linspace(0,1,10))
#     gl = GasLayer(10)
#     print(gl.get_arr())
#     print(c.foo(88))