import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import pickle
import os

if os.path.exists('Z:/users/wcapecch/'):
	direc = 'Z:/users/wcapecch/'
	proj_direc = 'Z:/PycharmProjects/LTXb-py/'
elif os.path.exists('//samba/wcapecch/'):
	direc = '//samba/wcapecch/'
	proj_direc = '//samba/wcapecch/PycharmProjects/LTXb-py/'


def field_line_pitch_not_constant_on_flux_surface():
	orbit = f'{proj_direc}neutral_beam/resonance/orbits/r10cm_rtan22cm_13kev_orbit.p'
	orbit_dat = pickle.load(open(orbit, 'rb'))
	out = orbit_dat['orb']['OUT']
	norm_psi = 0.8  # pick a flux contour to examine
	cs = plt.contour(out['X_XY'], out['Y_XY'], out['NORM_PSIXY'], levels=[norm_psi])
	for item in cs.collections:
		for i in item.get_paths():
			v = i.vertices  # vertices of contour
	xcont, ycont = v[:, 0], v[:, 1]
	xx, yy = out['X_XY'].flatten(), out['Y_XY'].flatten()
	
	# check our interpolation: zz below should be = norm_psi variable we used to define contour level
	zz = griddata((xx, yy), out['NORM_PSIXY'].flatten(), (xcont, ycont), method='cubic')
	
	bp = griddata((xx, yy), out['BP_XY'].flatten(), (xcont, ycont))
	bt = griddata((xx, yy), out['BPHI_XY'].flatten(), (xcont, ycont))
	line_pitch = bp / bt
	
	xmag = xcont - 0.4  # xcoords relative to mag axis
	theta = np.arctan2(ycont, xmag)
	
	plt.plot(out['RVES'], out['ZVES'], 'k')
	plt.plot(xcont, line_pitch)
	plt.xlabel('X (m)')
	plt.ylabel('Y (m)\nline_pitch')
	plt.show()


if __name__ == '__main__':
	field_line_pitch_not_constant_on_flux_surface()
