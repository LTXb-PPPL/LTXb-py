import matplotlib.pyplot as plt
import numpy as np
import matplotlib

clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
matplotlib.use('TkAgg')  # allows plotting in debug mode

direc = 'Z:/PycharmProjects/LTXb-py/npa/calibration/'
fp = f'{direc}test_data.txt'

with open(fp, 'r') as file:
	lines = file.readlines()
	dat = np.array([[float(a) for a in l.split('\t')] for l in lines])

for i in range(len(dat[0, :])):
	plt.plot(dat[:, i], label=f'sig{i}')
plt.legend()
plt.xlabel('"time"')
plt.ylabel('sig (V?)')
plt.show()

if __name__ == '__main__':
	pass
