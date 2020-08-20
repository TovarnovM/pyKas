# import os
# import sys
# from math import *
# wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
# wd = os.path.dirname(os.path.dirname(wd))
# sys.path.append(wd+"\\src")
# import numpy as np
# import pytest
# from pytest import approx

# from termodyn import DirectBallMany, get_optsmany_sample
# import matplotlib.pyplot as plt

# import pprint


# def _main1():
#     opts = get_optsmany_sample()
#     dbm = DirectBallMany(opts)
#     res = dbm.run()
#     plt.plot(res[:, 0], res[:, 1])
#     pprint.pprint(opts)
#     plt.show()

# if __name__ == "__main__":
#     _main1()
