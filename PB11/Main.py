from Invariants.Plot import *
from Pneumatics.PnInit import pn_create_layer
from OneVelocity.OvInit import ov_create_layer
from ElasticPiston.ElPistInit import el_pist_create_layer
# from Solver import solver
# from Pizdetc import solver
# from Test_gas_solver import solver
from Test_powd_solver import solver
# from Test_gas_powder_solver import solver
import time

# !!!!!!!!!!!!!!!!!!!!!!!! РАБОЧАЯ ПРОГРАММА !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
    # p_arr = [layers[-1].p[-1]]
    p_arr = layers[-1].p
    x_c_arr = layers[-1].x_c
    Vmax = 0
    pmax = 0
    results = []

    ind = 0

    while True:
        if (solver['geom'][-1][0] - layers[-1].x[-1]) <= 0.001:
            break
        if (Vmax - layers[-1].V[-1]) > 20:
            print('Замедляется')
            break
        if layers[0].time >= 0.1:
            print('Истекло время')
            break
        if pmax > 2e+9:
            print('Превысило давление')
            break
        # if (layers[-1].x[-1] - 5) >= 0.001:
        #     print('Вылетел из 15 метрового ствола')
        #     break

        layers1 = []
        tau_arr = [l.time_step() for l in layers]
        tau = solver['courant_number'] * min(tau_arr) # Вычисление шага по времени

        for i in range(len(layers)):
            if len(layers) == 1:
                layers1.append(layers[i].euler_step(layers, i, tau, flag='one'))
            elif i == 0:
                layers1.append(layers[i].euler_step(layers, i, tau, flag='left'))
                # layers[i] = layers1[i]
            elif i != 0 and i != (len(layers)-1):
                layers1.append(layers[i].euler_step(layers, i, tau, flag='inter'))
                # layers[i] = layers1[i]
            elif i == (len(layers)-1):
                if ind % 500 == 0:
                    layers1.append(layers[i].euler_step(layers, i, tau, flag='right'))
                else:
                    layers1.append(layers[i].euler_step(layers, i, tau, flag='right'))

        results.append(layers)
        layers = layers1

        time_arr.append(layers[-1].time)         # Добавление текущего шага по времени в список для графика
        V_arr.append(layers[-1].V[-1])           # Добавление текущей скорости поршня в список для графика
        x_arr.append(layers[-1].x[-1])           # Добавление текущей координаты поршня в список для графика
        # p_arr.append(layers[-1].p[-1])

        x_c_arr = np.vstack([x_c_arr, layers[-1].x_c])
        p_arr = np.vstack([p_arr, layers[-1].p])

        if layers[-1].V[-1] > Vmax:
            Vmax = layers[-1].V[-1]
        pmax = max([max(layer.p) for layer in layers])

        ind += 1
        if ind % 500 == 0:
            ro, u, e, z = layers[0].get_param(layers[0].q)
            print(layers[-1].time, Vmax, layers[0].V[-1], u[-1])
            # print(layers[1].p[-1])

    print("--- %s seconds ---" % (time.time() - start_time))

    print('Время вылета:', time_arr[-1], 'с')
    print('Скорость вылета:', V_arr[-1], 'м/с')
    print('Скорость вылета2:', layers[-1].u[-1], 'м/с')

    plot_one(time_arr, V_arr, 'График скорости снаряда от времени', 'Время', 'Скорость')
    plot_one(x_arr, V_arr, 'График скорости снаряда от времени', 'Координата', 'Скорость')
    plot_one(time_arr, p_arr, 'График давления на дно снаряда от времени', 'Время', 'Давление')
    plot_contours(x_c_arr, x_c_arr[0], time_arr, p_arr, 'Распределение давления', 'Координата', 'Время')
    plot_3d(x_c_arr, x_c_arr[0], time_arr, p_arr, 'Распределение давления на дно снаряда', 'Координата', 'Время')
    # plot_muzzle(results[0])


if __name__ == '__main__':
    main()
