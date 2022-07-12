import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav
import pickle
import os
import glob


def convert_sav_to_dict(fn):
	dat = readsav(fn)
	orb = {name: dat['orb'][0][name] for name in dat['orb'].dtype.names}
	for subrec in ['CIQ', 'OUT', 'AMB']:
		orb[subrec] = {name: orb[subrec][0][name] for name in orb[subrec].dtype.names}
	indict = {name: dat['in'][0][name] for name in dat['in'].dtype.names}
	orbit_dat = {'orb': orb, 'in': indict}
	pickle.dump(orbit_dat, open(f'{fn[:-4]}.p', 'wb'))
	print(f'pickled orbit: {fn}')


def test_pickled_orbit(fn):
	orbit_dat = pickle.load(open(fn, 'rb'))
	plt.plot(orbit_dat['orb']['R'], orbit_dat['orb']['Z'])
	plt.plot(orbit_dat['orb']['OUT']['RLIMIT'], orbit_dat['orb']['OUT']['ZLIMIT'])
	plt.plot(orbit_dat['orb']['OUT']['RVES'], orbit_dat['orb']['OUT']['ZVES'])
	plt.show()


if __name__ == '__main__':
	if os.path.exists('Z:/users/wcapecch/'):
		direc = 'Z:/users/wcapecch/'
		proj_direc = 'Z:/PycharmProjects/LTXb-py/'
	elif os.path.exists('//samba/wcapecch/'):
		direc = '//samba/wcapecch/'
		proj_direc = '//samba/wcapecch/PycharmProjects/LTXb-py/'
	pfns = glob.glob(f'{proj_direc}neutral_beam/resonance/orbits/*.p')
	test_pickled_orbit(pfns[0])

# orbit_dir = f'{proj_direc}neutral_beam/resonance/orbits/'
# orbits = glob.glob(orbit_dir + '*orbit.sav')
# convert_sav_to_dict(orbits[0])
