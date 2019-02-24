from ipyparallel import Client
from subprocess import Popen
import subprocess

_ipcontroller_process = None

class ClystyRx(object):
    def __init__(self, lbv):
        self.lbv = lbv



def is_controller_on():
    try:
        c = Client()
        return True
    except Exception as e:
        if str(e).startswith('Connection file'):
            return False
        print(e)
    raise Exception('А хз')

def run_controller():
    # global _ipcontroller_process
    # if _ipcontroller_process is None:
    #     _ipcontroller_process = Popen(['ipcontroller', "--ip='*'"], shell=True, stdin=subprocess.PIPE)
    global _ipcontroller_process
    _ipcontroller_process = Popen(['ipcontroller', "--ip='*'"], shell=True, stdin=subprocess.PIPE)
    return _ipcontroller_process
    # return Popen(['ipcontroller', "--ip='*'"])

def stop_controller():
    # if is_controller_on():
    #     pass
    _ipcontroller_process.communicate(subprocess)

def _forget_main():
    run_controller()
    if not is_controller_on():
        run_controller()
    print(is_controller_on())
    stop_controller()
    input()
    print(is_controller_on())

if __name__ == "__main__":
    pass