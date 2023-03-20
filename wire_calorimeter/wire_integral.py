import matplotlib.pyplot as plt
import numpy as np
from toolbox.helpful_stuff import make_patch_spines_invisible
from wire_calorimeter.wire_controls import WireControls

'''
given beam footprint (assumed gaussian), compute power intercepted by wire at various positions
todo: generalize to independent x, y gaussian
This is to model predicted dT for given wire dimensions
'''
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']

c = WireControls()
do_convergence_check = False


def xygauss_amp(x, y):
	return c.A * np.exp(-(y) ** 2 / (2. * c.sigy ** 2)) * np.exp(-(x) ** 2 / (2. * c.sigx ** 2))


def gauss_integral_check():
	xchk, ychk = np.linspace(-1000, 1000, endpoint=True, num=1000), np.linspace(-1000, 1000, endpoint=True, num=1000)
	dx, dy = xchk[1] - xchk[0], ychk[1] - ychk[0]
	ptot = 0.
	for x in xchk:
		for y in ychk:
			ptot += xygauss_amp(x, y) * dx * dy
	print(f'Ptot = {ptot * 1.e-3:.2f} (kW)')


# gauss_integral_check()  # OK


def ygauss_amp(x):
	# return amplitude of guassian in y given location x
	return c.A * np.exp(-(x) ** 2 / (2. * c.sigx ** 2))


def xgauss_amp(y):
	return c.A * np.exp(-(y) ** 2 / (2. * c.sigy ** 2))


def ygauss_int(x):
	# return indefinite integral of gaussian in y at location x
	return c.A * np.exp(-x ** 2 / (2. * c.sigx ** 2)) * np.sqrt(2. * np.pi * c.sigy ** 2)


if do_convergence_check:
	xw = -2.5  # x loc of wire
	fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
	nxarr = np.arange(10, 210, 10)  # number of bins across the wire
	amp = ygauss_amp(xw)
	power_computed = []
	gauss_integral = []
	for nx in nxarr:
		xx = np.linspace(xw - c.w_diam / 2., xw + c.w_diam / 2., nx + 1)  # make nx bins between bounds of wire
		dx, xmid = xx[1] - xx[0], [(xx[i + 1] + xx[i]) / 2. for i in np.arange(nx)]
		power2wire = [ygauss_int(xmid[i]) * dx for i in np.arange(nx)]
		power_computed.append(np.sum(power2wire))
		yy = np.linspace(-c.ygap / 2., c.ygap / 2., endpoint=True, num=nx + 1)
		dy, ymid = yy[1] - yy[0], np.array([(yy[i + 1] + yy[i]) / 2. for i in np.arange(nx)])
		z = amp * np.exp(-ymid ** 2 / (2 * c.sigy ** 2))
		gauss_integral.append(np.sum(z) * dy)
	ax1.plot(nxarr, power_computed, 'o-')
	ax2.plot(nxarr, gauss_integral, 'o-')
	ax1.axvline(c.nwire_bins, ls='--')
	ax2.axvline(c.ngauss_bins, ls='--')
	ax2.set_xlabel('nx')
	ax1.set_ylabel('power to wire [W]')
	ax2.set_ylabel('definite gauss integral')
	ax1.set_title('convergence check')
	plt.show()


# gauss_integral_check2()  # OK

def compute_fraction_to_wire_grid(wire_pos_arr=None, plot_changes_to_wires=False):
	# assume convergence is good, can check above
	if wire_pos_arr is None:
		xwire_pos_arr = (np.arange(c.nxwires) - (c.nxwires - 1) / 2) / c.nxwires * c.xgap  # equal gaps, 1/2 gap at edge
		ywire_pos_arr = (np.arange(c.nywires) - (c.nywires - 1) / 2) / c.nywires * c.ygap
	else:
		xwire_pos_arr, ywire_pos_arr = wire_pos_arr, wire_pos_arr
	xx = np.linspace(-c.xgap / 2., c.xgap / 2., endpoint=True, num=c.ngauss_bins + 1)  # bins along wire
	yy = np.linspace(-c.ygap / 2., c.ygap / 2., endpoint=True, num=c.ngauss_bins + 1)
	dx, xmid = xx[1] - xx[0], np.array([(xx[i + 1] + xx[i]) / 2. for i in np.arange(c.ngauss_bins)])
	dy, ymid = yy[1] - yy[0], np.array([(yy[i + 1] + yy[i]) / 2. for i in np.arange(c.ngauss_bins)])
	def_xgauss = np.sum(
		np.exp(-xmid ** 2 / (2 * c.sigx ** 2)) * dx)  # for horizontal wires; amp = 1, apply real amplitude later
	def_ygauss = np.sum(np.exp(-ymid ** 2 / (2 * c.sigy ** 2)) * dy)  # for vertical wires (unshadowed value)
	
	def gauss_integral_check2():
		nchk = 777
		ww = np.linspace(-c.xgap / 2., c.xgap / 2., endpoint=True, num=nchk)
		dw, wmid = ww[1] - ww[0], [(ww[i + 1] + ww[i]) / 2. for i in np.arange(nchk - 1)]
		power2wire = [ygauss_amp(wmid[i]) * def_ygauss * dw for i in np.arange(nchk - 1)]  # integrated_power*area
		power = np.sum(power2wire)
		print(f'Ptot = {power * 1.e-3:.2f} (kW)')
	
	vert_wire_power_computed = np.array([])
	horz_wire_power_computed, power_shadowed = np.array([]), np.array([])
	for wire_pos in xwire_pos_arr:  # compute power to vertical wires with NO shadowing
		ww = np.linspace(wire_pos - c.w_diam / 2., wire_pos + c.w_diam / 2., c.nwire_bins + 1)  # make nx bins btwn bounds of wire
		dw, wmid = ww[1] - ww[0], [(ww[i + 1] + ww[i]) / 2. for i in np.arange(c.nwire_bins)]
		power2wire = [ygauss_amp(wmid[i]) * def_ygauss * dw for i in np.arange(c.nwire_bins)]  # integrated_power*area
		vert_wire_power_computed = np.append(vert_wire_power_computed, np.sum(power2wire))
	for wire_pos in ywire_pos_arr:  # compute power to horizontal wires
		ww = np.linspace(wire_pos - c.w_diam / 2., wire_pos + c.w_diam / 2., c.nwire_bins + 1)  # make nx bins btwn bounds of wire
		dw, wmid = ww[1] - ww[0], [(ww[i + 1] + ww[i]) / 2. for i in np.arange(c.nwire_bins)]
		power2wire = [xgauss_amp(wmid[i]) * def_xgauss * dw for i in np.arange(c.nwire_bins)]  # integrated_power*area
		horz_wire_power_computed = np.append(horz_wire_power_computed, np.sum(power2wire))
	for wire_pos in xwire_pos_arr:  # compute shadowing per vertical wire
		pshadow = 0.  # shadowed power
		for ypos in ywire_pos_arr:
			pshadow += xygauss_amp(wire_pos, ypos) * c.w_diam ** 2  # assumes small wire diam so guass is ~const
		power_shadowed = np.append(power_shadowed, pshadow)
	vert_wire_power_corrected = vert_wire_power_computed - power_shadowed
	
	# compute predicted temperature change
	# P/A = rho*l*cp*dT/tpulse from Titus' paper
	# Q = m*c*dT
	ywire_mass = np.pi * (c.w_diam / 2.) ** 2 * c.ygap * c.rho  # [g]
	energy_to_ywires = vert_wire_power_corrected * c.tbeam  # [J]
	dtempy = energy_to_ywires / ywire_mass / c.c_sp  # [K]
	energy_to_ywires_uncorrected = vert_wire_power_computed * c.tbeam  # [J]
	dtempy_uncorrected = energy_to_ywires_uncorrected / ywire_mass / c.c_sp  # [K]
	xwire_mass = np.pi * (c.w_diam / 2.) ** 2 * c.xgap * c.rho
	energy_to_xwires = horz_wire_power_computed * c.tbeam
	dtempx = energy_to_xwires / xwire_mass / c.c_sp
	
	# compute predicted max temp range (center of center wire)
	# consider wire of length equal to its diameter (square area) at center of beam: power = A
	energy_to_center = c.A * c.tbeam * c.w_diam * c.w_diam
	mass_center = np.pi * (c.w_diam / 2.) ** 2 * c.w_diam * c.rho
	dtemp_center = energy_to_center / mass_center / c.c_sp  # [K]
	
	# compute predicted lengthening of wires due to increase temp
	dlx = c.xgap * dtempx * c.alpha  # [cm]
	dly = c.ygap * dtempy * c.alpha  # [cm]
	dly_uncorr = c.ygap * dtempy_uncorrected * c.alpha  # [cm]
	
	wire_intercept_fraction = np.sum([np.sum(horz_wire_power_computed), np.sum(vert_wire_power_corrected)]) / c.pbeam
	
	if plot_changes_to_wires:
		fs = 10
		fig, ax = plt.subplots()
		axr = ax.twinx()
		axrr = ax.twinx()
		axrr.spines["right"].set_position(("axes", 1.2))
		make_patch_spines_invisible(axrr)
		axrr.spines["right"].set_visible(True)
		axrr.yaxis.set_label_position('right')
		axrr.yaxis.set_ticks_position('right')
		axrr.set_ylabel('wire lengthening [mm]', fontsize=fs)
		ax.plot(ywire_pos_arr, horz_wire_power_computed * 1.e-3, 's-', label='horiz. wires', c=clrs[0])
		ax.plot(xwire_pos_arr, vert_wire_power_computed * 1.e-3, 'o-', label='vert. wires unshadowed', c=clrs[1])
		ax.plot(xwire_pos_arr, vert_wire_power_corrected * 1.e-3, 'o--', label='vert. wires shadowed', c=clrs[2])
		axr.plot(ywire_pos_arr, dtempx, 's-', c=clrs[0])
		axr.plot(xwire_pos_arr, dtempy_uncorrected, 'o-', c=clrs[1])
		axr.plot(xwire_pos_arr, dtempy, 'o--', c=clrs[2])
		axrr.plot(ywire_pos_arr, dlx * 10., 's-', c=clrs[0])  # 10x to convert cm to mm
		axrr.plot(xwire_pos_arr, dly_uncorr * 10., 'o-', c=clrs[1])
		axrr.plot(xwire_pos_arr, dly * 10., 'o-', c=clrs[2])
		ax.legend()
		ax.set_xlabel('wire position (cm)')
		ax.set_ylabel('power to wire (kW)')
		ax.set_title(f'Beam power to wires: {wire_intercept_fraction * 100:.2f}%, max dT: {dtemp_center:.1f} K')
		axr.set_ylabel('wire <dT> (K)')
		plt.tight_layout()
		plt.show()
	return wire_intercept_fraction

if __name__ == '__main__':
	pass
