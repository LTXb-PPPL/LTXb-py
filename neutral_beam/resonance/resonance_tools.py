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
	for k in ['R', 'PHI', 'Z', 'X', 'Y', 'VR', 'VPHI', 'VZ', 'BR', 'BPHI',
	          'BZ']:  # these have shape (20000,1) for some reason
		orb[k] = orb[k].flatten()
	out = orb['OUT']
	iz0 = np.where(abs(out['Y']) == min(abs(out['Y'])))[0][0]
	normpsi1d = out['NORM_PSIXY'][iz0, :]
	normpsi = out['NORM_PSIXY']
	if min(normpsi1d) == normpsi1d[0] or min(normpsi1d) == normpsi1d[-1]:  # Need to flip, not minimal on axis
		normpsi1d = -normpsi1d
		normpsi = -normpsi
		print(f'CAUTION: flipping PSI to put minimum on axis- might want to check this. File : {fn}')
	normpsi1d = normpsi1d - min(normpsi.flatten())  # remove any offset
	normpsi = normpsi - min(normpsi.flatten())
	psirmin = np.interp(min(out['RLIMIT']), out['X'], normpsi1d)
	psirmax = np.interp(max(out['RLIMIT']), out['X'], normpsi1d)
	assert (psirmin - psirmax) / np.mean(psirmin + psirmax) < 1.e-3  # psirmin should be = psirmax
	normpsi1d = normpsi1d / np.mean([psirmin, psirmax])  # now have psi= 0 -> 1 from axis to LCFS
	normpsi = normpsi / np.mean([psirmin, psirmax])
	out['NORM_PSI'] = normpsi1d
	out['NORM_PSIXY'] = normpsi
	orb['OUT'] = out
	orbit_dat = {'orb': orb, 'in': indict}
	pickle.dump(orbit_dat, open(f'{fn[:-4]}.pkl', 'wb'))
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
	pfns = glob.glob(f'{direc}PycharmProjects/data//orbits/*.pkl')
	test_pickled_orbit(pfns[0])

# orbit_dir = f'{direc}PycharmProjects/data/orbits/'
# orbits = glob.glob(orbit_dir + '*orbit.sav')
# convert_sav_to_dict(orbits[0])
