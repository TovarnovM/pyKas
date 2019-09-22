# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False

from gaslayer cimport GasLayer

cdef class GasPhase(GasLayer):
    pass

class MultiPhaseLayer(object):
    def __init__(self, layers):
        self.layers = layers

    def copy(self):
        res = MultiPhaseLayer({name:lr.copy() for name, lr in self.layers.items()})
        return res

    def  
    
    