import matplotlib.pyplot as plt
import numpy as np
from toolbox.particle_trajectory_tools import get_generic_path
from toolbox.helpful_stuff import ltx_above, ltx_limiter

clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']

xi, xo, yi, yo = ltx_above()
rlimiter, zlimiter, rminor_lim, theta_lim = ltx_limiter()
fig, (ax1, ax2) = plt.subplots(ncols=2)
ax1.plot(rlimiter, zlimiter)
ax2.plot(xi, yi)
ax2.plot(xo, yo)

for i in range(100):
	xp, yp, rp, zp, dp = get_generic_path(.4, 0., 0.)
	ax1.plot(rp, zp)
	ax2.plot(xp, yp)

plt.show()

if __name__ == '__main__':
	pass
