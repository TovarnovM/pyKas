#python setup_all.py build_ext --inplace
import glob
import os
import sys
import shutil
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy
from clean_stuff import cleanstuff


if __name__ == '__main__':
    true_mod_name='all'
    setup(
        name=true_mod_name,
        ext_modules = cythonize(["*.pyx"], annotate=True),
        include_dirs=[numpy.get_include()]
    )
    cleanstuff()