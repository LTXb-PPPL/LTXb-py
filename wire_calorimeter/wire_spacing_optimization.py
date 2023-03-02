import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

"""
This routine is aimed at answering the question of how to optimally space wires in an array to both maximize resolution
of gaussian distribution and minimize beam obstruction. Or rather maximize resolution for a given beam obstruction.
Two problems come to mind:
1) How do we define "resolution of a gaussian"
2) How do we systematically test different wire spacings?
thoughts on 1) are that given some noise level in our measurements, we can relate the signal to noise ratio generically
based simply on location from gaussian center. Then given our measurement locations and associated "uncertainties" based
on this SNR, we can determine some sort of likelihood that a fit to this data is representative (p-value?) of the real
distribution. Perhaps this p-value can act as our "good resolution" metric
thoughts on 2)- see mapping below: uses linear change to gap- could investigate quadratic or other
"""

clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']

num_wires = 7
half_gap = 12  # [cm] the aperture half-width
fig, ax = plt.subplots()
num_arrangements = 13
cumsum = np.sum(np.arange(1, num_wires))  # sum from 0:num_wires-1
ylim = half_gap * (num_wires - 1) / (num_wires ** 2 - num_wires - cumsum)  # when yo > ylim curve overshoots
yo_arr = np.linspace(0, ylim, num=num_arrangements + 1, endpoint=True)
yo_arr = yo_arr[1:-1]  # cut off 0 and ylim, others are equally spaced
m_arr = (half_gap - num_wires * yo_arr) / cumsum
for i in np.arange(len(yo_arr)):
	cs = np.append([0], np.cumsum(yo_arr[i] + m_arr[i] * np.arange(0, num_wires)))
	actual_wire_positions = cs[:-1]  # last point is at ygap, meant to ignore, included to make plot sensible
	ax.plot(np.arange(0, num_wires + 1), cs, 'o-', label=f'{m_arr[i]:.2f}')
ax.axhline(half_gap, c='k', ls='--')
ax.set_xlabel('wire #')
ax.set_ylabel('wire pos [cm]')
ax.legend(title='metric')

'''
we have generation method above to get actual_wire_position arrays for various configurations
need to compute a p-value or similar score for doing a "good job" at resolving gaussian
'''







plt.show()

if __name__ == '__main__':
	pass
