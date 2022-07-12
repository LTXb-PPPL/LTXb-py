import matplotlib.pyplot as plt
import numpy as np
import os
from helpful_stuff import read_eqdsk, ltx_limiter
from scipy.io import readsav

if os.path.exists('Z:/users/wcapecch/'):
	direc = 'Z:/users/wcapecch/'
	proj_direc = 'Z:/PycharmProjects/LTXb-py/'
elif os.path.exists('//samba/wcapecch/'):
	direc = '//samba/wcapecch/'
	proj_direc = '//samba/wcapecch/PycharmProjects/LTXb-py/'


def eqdsk_test():
	eqfn = f'{direc}/datasets/LTX_100981_468-1_75.eqdsk'
	out = read_eqdsk(eqfn)
	plt.contourf(out['x_xy'], out['y_xy'], out['psixy'], levels=25)
	plt.plot(out['rlimit'], out['zlimit'], 'k--')  # LCFS
	rlim, zlim, _, _ = ltx_limiter()
	plt.plot(rlim, zlim, 'k-')  # shell limiter
	plt.show()


def orbit_restore():
	orb_fn = f'{proj_direc}neutral_beam/resonance/idl_files/r10cm_rtan22cm_13kev_orbit.sav'
	orb = readsav(orb_fn)
	oorb = orb['orb']
	plt.plot(oorb.r[0], oorb.z[0])
	plt.plot(oorb.out[0].rlimit[0], oorb.out[0].zlimit[0], label='lim')
	plt.plot(oorb.out[0].rves[0], oorb.out[0].zves[0], label='ves')
	plt.legend()
	plt.show()
	a = 1


if __name__ == '__main__':
	# eqdsk_test()
	orbit_restore()
	pass
