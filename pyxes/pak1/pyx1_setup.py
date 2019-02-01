#python pyx1_setup.py build_ext --inplace
import glob
import os
import shutil
from distutils.core import setup
from Cython.Build import cythonize

if __name__ == '__main__':
    true_mod_name='pyx1'
    setup(
        name='pyx1',
        ext_modules = cythonize("pyx1.pyx", language="c++")
    )

    try:
        f2 = glob.glob(true_mod_name+'.pyd')
        os.remove(f2[0])
    except Exception as e:
        print(e)

    f = glob.glob(true_mod_name+'*.pyd')
    os.rename(f[0], true_mod_name+'.pyd')
    # os.remove(pyx_name)
    os.remove(true_mod_name+'.cpp')
    shutil.rmtree('build')