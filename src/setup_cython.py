#python setup_cython.py build_ext --inplace
import glob
import os
import sys
import shutil
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from clean_stuff import cleanstuff
import numpy


if __name__ == '__main__':
    true_mod_name='all'
    extentions = [
        Extension("tube", ["include\\tube\\tube.pyx"] ,include_dirs=[numpy.get_include()]),
        Extension("gaslayer", ["include\\gaslayer\\gaslayer.pyx"], library_dirs=['\\include\\tube\\'] ),
    ]
    ext_modules = cythonize(extentions, annotate=True)
    setup(
        name=true_mod_name,
        ext_modules=ext_modules,
    )
    cleanstuff()