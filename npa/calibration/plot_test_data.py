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
# fp = f'{hdirec}npa/calibration/anpa_calib_test.dat'
# fp = f'{hdirec}npa/calibration/first_scan.dat'
# fp = f'{hdirec}npa/calibration/second_scan.dat'
fp = f'{hdirec}npa/calibration/fourth_scan.dat'

with open(fp, 'r') as file:
	lines = file.readlines()
	dat = np.array([[float(a) for a in l.split('\t')] for l in lines])

if len(dat[0, :]) > len(dat[:, 0]):
	dat = dat.transpose()

# vsrc = dat[:, 14]
# [02,03,04,05,06,07,08,09,10,11,12,13,14]  # AI
# [04,05,00,06,07,00,00,00,00,08,09,10,rr]  # T
# [00,01,02,03,04,05,06,07,08,09,10,11,12]  # dat index
fig, ax = plt.subplots(nrows=2, ncols=4, sharex='col', sharey='row')
ax[0, 0].plot(dat[:, 0], label='T4')
ax[0, 1].plot(dat[:, 1], label='T5')
ax[0, 2].plot(dat[:, 3], label='T6')
ax[0, 3].plot(dat[:, 4], label='T7')
ax[1, 0].plot(dat[:, 9], label='T8')
ax[1, 1].plot(dat[:, 10], label='T9')
ax[1, 2].plot(dat[:, 11], label='T10')
ax[1, 3].plot(dat[:, 12], label='ross')

for aa in ax:
	for a in aa:
		a.legend()
plt.xlabel('"time"')
plt.ylabel('sig (V)')
plt.show()

if __name__ == '__main__':
	pass
