from Invariants.Plot import plot_one
from Pneumatics.PnInit import pn_create_layer
from OneVelocity.OvInit import ov_create_layer
from Tests.Dict import solver


def main():
    layers = []
    for i in range(len(solver['grids'])):
        if solver['grids'][i]['type'] == 'gas':
            layers.append(pn_create_layer(i))
        elif solver['grids'][i]['type'] == 'powder':
            layers.append(ov_create_layer(i))

    time_arr = [0]              # Список времени для графика
    V_arr = [0]                 # Список скорости для графика
    x_arr = [layers[-1].x[-1]]           # Список координаты для графика

    while (solver['borders'][-1]['x'] - layers[-1].x[-1]) >= 0.001:
        tau_arr = [l.time_step() for l in layers]
        tau = solver['courant_number'] * min(tau_arr) # Вычисление шага по времени
        layer1 = layer.euler_step(layer, tau, p_right=p_atm)    # Шаг по времени
        layer = layer1
        time_arr.append(layer.time)         # Добавление текущего шага по времени в список для графика
        V_arr.append(layer.V[-1])           # Добавление текущей скорости поршня в список для графика
        x_arr.append(layer.x[-1])           # Добавление текущей координаты поршня в список для графика

    print('Время вылета:', time_arr[-1], 'с')
    print('Скорость вылета:', V_arr[-1], 'м/с')

    plot_one(time_arr, V_arr, 'График скорости поршня от времени', 'Время', 'Скорость')
    plot_one(x_arr, V_arr, 'График скорости поршня от координаты', 'Координата', 'Скорость')
    plot_one(time_arr, x_arr, 'График координаты поршня от времени', 'Время', 'Координата')


if __name__ == '__main__':
    main()
