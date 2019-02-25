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
from shutil import copyfile


def compile_bicycle(packages):
    true_mod_name='all'
    incl_dir = "include"
    wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
    wd = os.path.dirname(wd)

    pxds = []
    for pn in packages:
        depends = packages[pn].get("depends")
        if depends:
            for dep in depends:
                pxd = f'{wd}\\{incl_dir}\\{dep}\\{dep}.pxd'
                print(f'Finding pxd ... {pxd}...', end='')
                if os.path.isfile(pxd):
                    print('FIND')
                    todir = f'{wd}\\{incl_dir}\\{pn}'
                    print(f'Finding to dir for {dep}.pxd...{todir}...', end='')
                    if os.path.isdir(todir):
                        print('FIND')
                        pxds.append((pxd, f'{todir}\\{dep}.pxd'))
                    else:
                        raise AttributeError(f'Папка {todir} не найженf')
                else:
                    raise AttributeError(f'Файл {pxd} не найжен')
    try:
        for src, dst in pxds:
            copyfile(src, dst)
        
        extentions = [
            Extension( 
                pn, 
                [f'{wd}\\{incl_dir}\\{pn}\\{pn}.pyx'], 
                **packages[pn]["Extention_kwargs"])
            for pn in packages]

        ext_modules = cythonize(extentions, annotate=True)
        setup(
            name=true_mod_name,
            ext_modules=ext_modules,
        )

    finally:
        for dont_del, must_del in pxds:
            try:
                print(f'removing {must_del}...',end='')
                os.remove(must_del)
                print('Done')
            except Exception as e:
                print('! ERROR !')
                print(e)
        


if __name__ == '__main__':
    packages = {
            "tube": {
                "Extention_kwargs":  { "include_dirs":[numpy.get_include()] }
            }, 
            
            "gaslayer": {
                "Extention_kwargs":{},
                "depends": ["tube"]
            }    
        }

    cleanstuff()
    compile_bicycle(packages)
    cleanstuff()
        
        # true_mod_name='all'
        # extentions = [
        #     Extension("tube", ["include\\tube\\tube.pyx"] ,include_dirs=[numpy.get_include()]),
        #     Extension("gaslayer", ["include\\gaslayer\\gaslayer.pyx"], library_dirs=['\\include\\tube'] ),
        # ]
        # ext_modules = cythonize(extentions, annotate=True)
        # setup(
        #     name=true_mod_name,
        #     ext_modules=ext_modules,
        #     library_dirs=['\\include\\tube']
        # )   