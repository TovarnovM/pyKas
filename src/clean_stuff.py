import os
import sys
import glob
import shutil
import pathlib
import itertools as it

def multiple_file_types(*patterns):
    return it.chain.from_iterable(glob.iglob(pattern) for pattern in patterns)

def cleanstuff():
    try:
        print("Try remove 'build\\' folder...", end=' ')
        shutil.rmtree('build')
        print("Done")
    except Exception as e:
        print(e)
    
    for file_cpp in multiple_file_types('*.cpp', "*.c"):
        try:
            print(f"Try remove '{file_cpp}' file...", end=' ')
            os.remove(file_cpp)
            print("Done")
        except Exception as e:
            print(e)

    for file_pyd in glob.glob('*.*.pyd'): 
        try:
            print(f"Try rename '{file_pyd}' file...", end=' ')
            names = file_pyd.split('.')
            good_name = names[0]+'.'+names[-1]
            if pathlib.Path(good_name).is_file():
                os.remove(good_name)
            os.rename(file_pyd, good_name)
            print("Done")
        except Exception as e:
            print(e)   

if __name__ == "__main__":
    cleanstuff()