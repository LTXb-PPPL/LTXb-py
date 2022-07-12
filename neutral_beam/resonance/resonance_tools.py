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
	a=1
	# pickle.dump(orbit_dat, open('', 'wb'))
	# print(f'pickled orbit: {fn}')

	
if __name__ == '__main__':
	if os.path.exists('Z:/users/wcapecch/'):
		direc = 'Z:/users/wcapecch/'
		proj_direc = 'Z:/PycharmProjects/LTXb-py/'
	elif os.path.exists('//samba/wcapecch/'):
		direc = '//samba/wcapecch/'
		proj_direc = '//samba/wcapecch/PycharmProjects/LTXb-py/'
	orbit_dir = f'{proj_direc}neutral_beam/resonance/orbits/'
	orbits = glob.glob(orbit_dir + '*orbit.sav')
	convert_sav_to_dict(orbits[0])