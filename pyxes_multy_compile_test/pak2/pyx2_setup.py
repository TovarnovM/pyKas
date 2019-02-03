#python pyx2_setup.py build_ext --inplace
import glob
import os
import sys
import shutil
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension


if __name__ == '__main__':
    sys.path.append(os.path.dirname(sys.path[0])+"\\pak1")
    true_mod_name='pyx2'
    setup(
        name=true_mod_name,
        ext_modules = cythonize(["..\\pak1\\pyx1.pyx", true_mod_name+".pyx"], language="c++")
    )

    # try:
    #     f2 = glob.glob(true_mod_name+'.pyd')
    #     os.remove(f2[0])
    # except Exception as e:
    #     print(e)

    # f = glob.glob(true_mod_name+'*.pyd')
    # os.rename(f[0], true_mod_name+'.pyd')
    # # os.remove(pyx_name)
    # os.remove(true_mod_name+'.cpp')
    # shutil.rmtree('build')