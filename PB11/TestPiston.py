from Invariants.Plot import *
from Pneumatics.PnInit import pn_create_layer
from OneVelocity.OvInit import ov_create_layer
from ElasticPiston.ElPistInit import el_pist_create_layer
from Solver import solver
import time


def main():
    start_time = time.time()
    layers = []
    for i in range(len(solver['grids'])):
        if solver['grids'][i]['type'] == 'gas':
            layers.append(pn_create_layer(i, solver))
        elif solver['grids'][i]['type'] == 'powder':
            layers.append(ov_create_layer(i, solver))
        elif solver['grids'][i]['type'] == 'piston':
            layers.append(el_pist_create_layer(i, solver))

    time_arr = [0]              # Список времени для графика
    V_arr = [0]                 # Список скорости для графика
    x_arr = [layers[-1].x[-1]]           # Список координаты для графика
    p_arr = [layers[-1].p[-1]]

    results = []

    while (solver['geom'][-1][0] - layers[-1].x[-1]) >= 0.001:
        layers1 = []
        tau_arr = [l.time_step() for l in layers]
        tau = solver['courant_number'] * min(tau_arr) # Вычисление шага по времени

        for i in range(len(layers)):
            if len(layers) == 1:
                layers1.append(layers[i].euler_step(layers[i], tau, p_left=4*10**7, p_right=101325))
            elif i == 0:
                layers1.append(layers[i].euler_step(layers[i], tau, p_right=layers[i+1].p[0]))
            elif i != 0 and i != (len(layers)-1):
                layers1.append(layers[i].euler_step(layers[i], tau, p_left=layers[i-1].p[-1], p_right=layers[i+1].p[0]))
            elif i == (len(layers)-1):
                layers1.append(layers[i].euler_step(layers[i], tau, p_left=layers[i-1].p[-1]))

        results.append(layers)
        layers = layers1

        time_arr.append(layers[-1].time)         # Добавление текущего шага по времени в список для графика
        V_arr.append(layers[-1].V[-1])           # Добавление текущей скорости поршня в список для графика
        x_arr.append(layers[-1].x[-1])           # Добавление текущей координаты поршня в список для графика
        p = max(layers[-1].p)
        p_arr.append(p)

    print("--- %s seconds ---" % (time.time() - start_time))

    print('Время вылета:', time_arr[-1], 'с')
    print('Скорость вылета:', V_arr[-1], 'м/с')

    plot_one(time_arr, V_arr, 'График скорости снаряда от времени', 'Время', 'Скорость')
    plot_one(time_arr, p_arr, 'График давления на дно снаряда от времени', 'Время', 'Давление')
    # plot_one(x_arr, V_arr, 'График скорости поршня от координаты', 'Координата', 'Скорость')
    # plot_one(time_arr, x_arr, 'График координаты поршня от времени', 'Время', 'Координата')

    # plot_contours(x_c_arr, x_c_arr[0], time_arr, p_arr, 'Распределение давления', 'Координата', 'Время')
    # plot_3d(x_c_arr, x_c_arr[0], time_arr, p_arr, 'Распределение давления', 'Координата', 'Время')

    plot_muzzle(results[0])


if __name__ == '__main__':
    main()
