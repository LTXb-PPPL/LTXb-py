import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from wire_calorimeter.wire_controls import WireControls
from wire_calorimeter.wire_integral import compute_fraction_to_wire_grid

# todo: REWORK wire_integral to take in wire positions (assume same positioning in x and y) and return relevant values (including intercepted fraction) then import and use here for various wire spacing configurations
# todo: Once we have intercept fraction for each spacing metric, compute "GOODNESS OF FIT" or something for each spacing metric

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

plot_wire_locations = True
plot_frac_v_wire_locations = True
if plot_wire_locations:
	fig, ax = plt.subplots()
if plot_frac_v_wire_locations:
	fig2, ax2 = plt.subplots()

c = WireControls()
cumsum = np.sum(np.arange(1, c.nxwires))  # sum from 0:num_wires-1
ylim = c.half_gap * (c.nxwires - 1) / (c.nxwires ** 2 - c.nxwires - cumsum)  # when yo > ylim curve overshoots
# yo_arr is array of gap between center and first wire
yo_arr = np.linspace(0, ylim, num=c.num_arrangements + 1, endpoint=True)
yo_arr = yo_arr[1:-1]  # cut off 0 and ylim, others are equally spaced
m_arr = (c.half_gap - c.nxwires * yo_arr) / cumsum  # slope of linear change in wire spacing: metric defining spacing


def wire_placements(metric):
	yo = -(metric * cumsum - c.half_gap) / c.nxwires
	cs = np.append([0], np.cumsum(yo + metric * np.arange(0, c.nxwires)))
	actual_wire_positions = cs[:-1]
	if plot_wire_locations:
		ax.plot(np.arange(0, c.nxwires + 1), cs, 'o-', label=f'{metric:.2f}')
	return actual_wire_positions


frac_arr = np.array([])
for m in m_arr:
	wp = wire_placements(m)  # one sided wire placements
	frac = compute_fraction_to_wire_grid(np.append(-wp[::-1], wp[1:]))
	frac_arr = np.append(frac_arr, frac)
if plot_frac_v_wire_locations:
	ax2.plot(m_arr, frac_arr*100., 'o-')
	ax2.set_xlabel('metric')
	ax2.set_ylabel('intercept (%)')

if plot_wire_locations:
	ax.axhline(c.half_gap, c='k', ls='--')
	ax.set_xlabel('wire #')
	ax.set_ylabel('wire pos [cm]')
	ax.legend(title='metric')

plt.show()

if __name__ == '__main__':
	pass
