import os
import sys
from math import *
wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
wd = os.path.dirname(os.path.dirname(wd))
sys.path.append(wd+"\\src\\")
sys.path.append(wd+"\\src\\include\\godunov")
from riman_python import plot_distrs, get_distrs_to_time
import numpy as np
import pytest
from pytest import approx

from tube import Tube
from gaslayer import GasLayer, GasEOS, GridStrecher, GasFluxCalculator
import time

def test_adapt():
    gs = GridStrecher()
    n = 100#np.random.randint(10,100)
    for i in range(1000):
        n = np.random.randint(10,1000)
        xx1 = np.random.uniform(-1000,2000)
        l = np.random.uniform(0.001, 2000)
        xs = np.linspace(xx1, xx1+l, n)
        xs_borders = np.linspace(xx1, xx1+l, n+1)

        ys = 100000*np.sin(xs)  
        xs_adapted = np.zeros_like(xs_borders)

        gs.adaptine_borders(xs_borders, ys, xs_adapted)
        assert xs_adapted[0] == approx(xs_borders[0])
        assert xs_adapted[-1] == approx(xs_borders[-1])
        for x1, x2 in zip(xs_adapted[:-1], xs_adapted[1:]):
            assert x2 > x1


def _plot():
    import matplotlib.pyplot as plt

    n = 1000
    xs = np.linspace(-10, 10, n)
    xs += np.random.uniform(-0.01,0.01,n)

    ys = np.sin(xs/2) + 0.3*np.random.normal(size=n)
    plt.plot(xs, ys)    
    
    ys_smoothed = np.zeros_like(ys)

    gs = GridStrecher()
    gs.smooth_arr(xs, ys, ys_smoothed,0.01)
    plt.plot(xs, ys_smoothed)
    
    gs.smooth_arr(xs, ys, ys_smoothed,0.1)
    plt.plot(xs, ys_smoothed)

    gs.smooth_arr(xs, ys, ys_smoothed,0.2)
    plt.plot(xs, ys_smoothed)

    plt.show()

def _plot2():
    import matplotlib.pyplot as plt

    n = 1000
    xs = np.linspace(-10, 10, n)
    xs_borders = np.linspace(-10, 10, n+1)

    ys = 100000*np.sin(xs)
    plt.plot(xs, ys)    
    
    xs_adapted = np.zeros_like(xs_borders)

    gs = GridStrecher()
    gs.adaptine_borders(xs_borders, ys, xs_adapted)
    plt.scatter(xs_adapted, np.zeros_like(xs_adapted))

    plt.show()

def _plot3():
    import matplotlib.pyplot as plt
    # из статьи ОДНОМЕРНЫЕ ЗАДАЧИ ГАЗОВОЙ ДИНАМИКИ И ИХ РЕШЕНИЕ  ПРИ ПОМОЩИ РАЗНОСТНЫХ СХЕМ ВЫСОКОЙ  РАЗРЕШАЮЩЕЙ СПОСОБНОСТИ 
    #           LEFT          RIGHT              
    #      ro       u       p       ro      u       p       t
    v = [( 1,       0,      1,      0.125,  0,      0.1,    0.15  ),
        #  ( 0.445,   0.698,  3.528,  0.5,    0,      0.571,  0.05  ),
        #  ( 1,       1,      1,      0.5,    0,      0.571,  0.05  )
         ]
    for ro_1, u_1, p_1, ro_2, u_2, p_2, t in v:
        d=  {
            'p_1' : p_1,          # давление слева
            'ro_1' : ro_1,         # плотность слева   
            'u_1' : u_1,          # скорость слева      
            'p_2' : p_2,        # давление справа     
            'ro_2': ro_2,      # плотность справа    
            'u_2' : u_2,          # скорость справа        
            'p_0' : 0,          # параметр в ур-ии состояния        
            'gamma' : 1.4,      # параметр в ур-ии состояния          
            'c_0' : 0,          # параметр в ур-ии состояния     
            'eps_F':1e-6,       # точность определения решения           
            'n_iter_max':100,   # максимальное количество итераций при определении решения              
            'x0' : 0.5,         # положение разрыва в момент t=0        
            'ts': [t],        # времена, в которых нужно построить графики распределения параметров         
            'n': 10000          # кол-во точек, в которых ищутся параметры волны разрежения         
        } 
    # plot_distrs(**d)
    # print(get_distrs_to_time(d['ts'][0], **d))

    x1 = 0
    x2 = 1
    def init_foo(x, *args):
        if x <= 0.5:
            return d['ro_1'], d['p_1'], d['u_1']
        return d['ro_2'], d['p_2'], d['u_2']
    
    gas_eos = GasEOS(gamma=1.4, kind=1)
    grid_strecher = GridStrecher(strech_type=1,st2_window_part=0.05, st2_adapt_prop=0.1)
    gas_flux_calculator = GasFluxCalculator()
    tube = Tube([0,1], [0.1, 0.1])
    n_cells = 300
    layer = GasLayer(n_cells, tube=tube, gasEOS=gas_eos, flux_calculator=gas_flux_calculator, grid_strecher=grid_strecher, n_qs=3)
    layer.xs_borders[:] = np.linspace(x1, x2, n_cells+1)
    layer.init_ropue_fromfoo(init_foo)


    # d = layer.to_dict()


    # grid_strecher.smooth_arr(layer.xs_cells, layer.es, layer.es, grid_strecher.st2_window_part)
    # grid_strecher.adaptine_borders(layer.xs_borders, layer.es, layer.xs_borders)
    # layer.init_ropue_fromfoo(init_foo)
    def plot_layer(lr):
        lr_d = layer.to_dict()
        plt.scatter(lr_d['xs_cells'], lr_d['ros'])
        plt.scatter(lr_d['xs_cells'], lr_d['ps'])
        plt.scatter(lr_d['xs_cells'], lr_d['us'])
        plt.scatter(lr_d['xs_cells'], lr_d['es'])
        # plt.scatter(lr_d['xs_borders'], lr_d['fluxes'][2])

        dd = get_distrs_to_time(layer.time, **d)
        plt.plot(dd['xs'], dd['ros'], label=r'$\rho$')
        plt.plot(dd['xs'], dd['ps'], label=r'$p$')
        plt.plot(dd['xs'], dd['us'], label=r'$u$')
        plt.plot(dd['xs'], dd['es'], label=r'$e$')
        plt.grid(True)
        plt.legend()
        plt.show()
    layer.init_taus_acustic()
    
    t0 = time.time()
    while layer.time < d['ts'][0]:
        # d['ts'][0] = 0.001
        # layer.time = 0.01
        # break
        tau = layer.get_tau_min()*0.5
        layer1 = layer.step_simple(tau, 0,0)
        if layer1 == layer:
            break
        # layer1 = layer.step_Godunov_corrector2(layer1, 0, 0)
        # if layer1 == layer:
        #     break
        layer = layer1
    print(time.time() - t0)  


    #(3.55189 + 3.40513 + 3.33072 + 3.44629 + 3.717145 + 3.93337 + 3.94355 + 3.706031+3.54699+3.4006)  
        # plot_layer(layer)
    print(layer.time)
    # print(layer.to_dict())
    plot_layer(layer)

if __name__ == "__main__":
    _plot3()
    # _plot()
    # _plot2()