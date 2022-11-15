import matplotlib.pyplot as plt
import numpy as np
import re
from toolbox.helpful_stuff import read_eqdsk2


def check_lim_file(fn, eqdsks=[]):
	f = open(fn, 'r')
	form1 = r'-?\d+.\d+E[+-]\d+'
	form2 = r'-?\d+.\d+e[+-]\d+'
	for i in np.arange(24):  # skip stuff in first 24 lines
		f.readline()
	nn = int(re.findall(r'\d+', f.readline())[0])
	nrows = np.ceil(nn / 6.)  # num rows to scan
	rlim, zlim = [], []
	for i in np.arange(nrows):
		if i == 0:  # do quick check for which form to use
			line = f.readline()
			if len(re.findall(form1, line)) > len(re.findall(form2, line)):
				form = form1
			else:
				form = form2
			rlim.extend([float(n) for n in re.findall(form, line)])  # already read line
		else:
			rlim.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		zlim.extend([float(n) for n in re.findall(form, f.readline())])
	f.close()
	
	for eq in eqdsks:
		dat = read_eqdsk2(eq)
		plt.plot(dat['rzout'][0, :], dat['rzout'][1, :])
	plt.plot(rlim, zlim, label=fn.split('/')[-1])
	plt.xlabel('R (m)')
	plt.ylabel('Z (m)')
	plt.tight_layout()


if __name__ == '__main__':
	# check_lim_file('Z:/transp/t105795/SAFE105795.LIM',
	#                eqdsks=[f'Y:/reconstructions/eqdsk/105795/run_{num}/Psitri.eqdsk' for num in np.arange(1, 24)])
	check_lim_file('Z:/transp/t105795/LTX105795.LIM')
	check_lim_file('Z:/transp/t105795/CAP105795.LIM')
	plt.legend()
	plt.show()
