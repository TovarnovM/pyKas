#python setup_multy.py build_ext --inplace
import glob
import os
import sys
import shutil
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension


if __name__ == '__main__':
    true_mod_name='pyx23'
    extentions = [
        Extension("pyx1", ["pak1\\pyx1.pyx"]),
        Extension("pyx2", ["pak2\\pyx2.pyx"], include_dirs=["pak1\\"])
    ]
    setup(
        name=true_mod_name,
        ext_modules = cythonize(extentions, language="c++", annotate=True)
    )