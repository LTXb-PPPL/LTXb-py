import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import os


def get_home_direc():
    pppl = 'Z:/PycharmProjects/LTXb-py/'
    home = 'C:/Users/willi/PycharmProjects/LTXb-py/'
    if os.path.exists(pppl):
        return pppl
    elif os.path.exists(home):
        return home
    else:
        print('no home direc found')
        return None


clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
matplotlib.use('TkAgg')  # allows plotting in debug mode

hdirec = get_home_direc()
# fp = f'{hdirec}npa/calibration/test_data.txt'
fp = f'{hdirec}npa/calibration/anpa_calib_test.dat'

with open(fp, 'r') as file:
    lines = file.readlines()
    dat = np.array([[float(a) for a in l.split('\t')] for l in lines])

if len(dat[0, :]) > len(dat[:, 0]):
    dat = dat.transpose()

for i in range(len(dat[0, :])):
    plt.plot(dat[:, i], label=f'sig{i}')
plt.legend()
plt.xlabel('"time"')
plt.ylabel('sig (V)')
plt.show()

if __name__ == '__main__':
    pass
