import matplotlib.pyplot as plt
import numpy as np
import random

clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']


def get_generic_path(x0, y0, z0):
	# This is ported from IDL routine /home/mstfit/mstfit/user_defined_routines/fast_ion/get_generic_div_path.pro
	# However this only takes 2 rotations and returns path at random orientation from starting point
	x1 = np.linspace(0., 1.5, num=1501)  # pt every mm
	dp = x1[1]-x1[0]
	phi, theta = random.random() * 2. * np.pi, random.random() * np.pi - np.pi/2.  # midplane, vertical angles
	A1 = np.arctan2(-y0, -x0) - phi  # Start with x1 in x^ direction, then rotate to machine CL then back by phi
	A2 = -theta  # positive theta angles up relative to midplane but requires negative rotation about y^
	xp = x0 + x1 * np.cos(A2) * np.cos(A1)
	yp = y0 + x1 * np.cos(A2) * np.sin(A1)
	zp = z0 - x1 * np.sin(A2)
	rp = np.sqrt(xp ** 2 + yp ** 2)
	return xp, yp, rp, zp, dp


if __name__ == '__main__':
	path = get_generic_path(.4, 0., 0.)
	pass
