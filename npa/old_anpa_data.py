import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

matplotlib.use('TkAgg')  # allows plotting in debug mode
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']


def get_old_anpa_data():
    # data from shot
    direc = 'C:/Users/willi/OneDrive/Desktop/work_stuff/NPA/npa_code/'
    dz = pd.read_csv(f'{direc}old_anpa_data_z.csv', header=None)
    dx = pd.read_csv(f'{direc}old_anpa_data_x.csv', header=None)
    dy = pd.read_csv(f'{direc}old_anpa_data_y.csv', header=None)
    x, y = dx.values.flatten(), dy.values.flatten()
    z = dz.values
    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex='col')
    ax1.contourf(x, y, z.transpose(), levels=np.linspace(0, .02))
    for i in range(11):
        ax2.plot(x, z[:, i], label=f'H{i + 1}')
    ax2.set_xlabel('time')

    ax1.set_ylim((0, 50))
    ax2.set_ylim((0, .05))
    a = 1
    return None


if __name__ == '__main__':
    dat = get_old_anpa_data()
