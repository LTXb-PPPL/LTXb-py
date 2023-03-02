import matplotlib.pyplot as plt
import numpy as np
from toolbox.helpful_stuff import make_patch_spines_invisible

'''
given beam footprint (assumed gaussian), compute power intercepted by wire at various positions
todo: generalize to independent x, y gaussian
This is to model predicted dT for given wire dimensions
'''
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']

# controls
nwire_bins = 100
ngauss_bins = 100
do_convergence_check = False

# wire layout
xgap, ygap = 24, 24  # [cm] aperture width in x, y (NOTE: 32 cm is TAE design)
nxwires, nywires = 7, 7  # num wires in each direction. xwires refers to x location of vertical wires

# beam characteristics
fwhmx = 8.  # cm
fwhmy = 9.  # cm
sigx = fwhmx / (2. * np.sqrt(2. * np.log(2)))  # [cm]
sigy = fwhmy / (2. * np.sqrt(2. * np.log(2)))  # [cm]
pbeam = 400.e3  # [W] total beam power
tbeam = 10.e-3  # [s] beam pulse duration

# wire characteristics: TAE design, 0.5mm diam wires 32 cm long
w_diam = .05  # [cm]

# Tungsten characteristics
c_mol = 24.27  # [J/mol/K] molar heat capacity
rho = 19.25  # [g/cm^3]
amu = 183.84  # [g/mol]
c_sp = c_mol / amu  # [J/g/K] specific heat capacity
k_con = 175.  # [W/m/K]  # thermal conductivity at 20C
alpha = 4.3e-6  # [K^-1] thermal expansion coeff at 20C

# 2d gaussian form- no covariance
# z = A * np.exp(-(x - xw) ** 2 / (2. * sigx ** 2) - (y - yw) ** 2 / (2. * sigy ** 2))
# so for given values of x, xw, sigx then z becomes 1d gaussian in y
# total integral: V=2*pi*A*sigx*sigy=total beam power
A = pbeam / (2. * np.pi * sigx * sigy)


def xygauss_amp(x, y):
	return A * np.exp(-(y) ** 2 / (2. * sigy ** 2)) * np.exp(-(x) ** 2 / (2. * sigx ** 2))


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
	return A * np.exp(-(x) ** 2 / (2. * sigx ** 2))


def xgauss_amp(y):
	return A * np.exp(-(y) ** 2 / (2. * sigy ** 2))


def ygauss_int(x):
	# return indefinite integral of gaussian in y at location x
	return A * np.exp(-x ** 2 / (2. * sigx ** 2)) * np.sqrt(2. * np.pi * sigy ** 2)


if do_convergence_check:
	xw = -2.5  # x loc of wire
	fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
	nxarr = np.arange(10, 210, 10)  # number of bins across the wire
	amp = ygauss_amp(xw)
	power_computed = []
	gauss_integral = []
	for nx in nxarr:
		xx = np.linspace(xw - w_diam / 2., xw + w_diam / 2., nx + 1)  # make nx bins between bounds of wire
		dx, xmid = xx[1] - xx[0], [(xx[i + 1] + xx[i]) / 2. for i in np.arange(nx)]
		power2wire = [ygauss_int(xmid[i]) * dx for i in np.arange(nx)]
		power_computed.append(np.sum(power2wire))
		yy = np.linspace(-ygap / 2., ygap / 2., endpoint=True, num=nx + 1)
		dy, ymid = yy[1] - yy[0], np.array([(yy[i + 1] + yy[i]) / 2. for i in np.arange(nx)])
		z = amp * np.exp(-ymid ** 2 / (2 * sigy ** 2))
		gauss_integral.append(np.sum(z) * dy)
	ax1.plot(nxarr, power_computed, 'o-')
	ax2.plot(nxarr, gauss_integral, 'o-')
	ax1.axvline(nwire_bins, ls='--')
	ax2.axvline(ngauss_bins, ls='--')
	ax2.set_xlabel('nx')
	ax1.set_ylabel('power to wire [W]')
	ax2.set_ylabel('definite gauss integral')
	ax1.set_title('convergence check')
	plt.show()

# assume convergence is good, can check above
xwire_pos_arr = (np.arange(nxwires) - (nxwires - 1) / 2) / nxwires * xgap  # equal gaps, 1/2 gap at edge
ywire_pos_arr = (np.arange(nywires) - (nywires - 1) / 2) / nywires * ygap
xx = np.linspace(-xgap / 2., xgap / 2., endpoint=True, num=ngauss_bins + 1)  # bins along wire
yy = np.linspace(-ygap / 2., ygap / 2., endpoint=True, num=ngauss_bins + 1)
dx, xmid = xx[1] - xx[0], np.array([(xx[i + 1] + xx[i]) / 2. for i in np.arange(ngauss_bins)])
dy, ymid = yy[1] - yy[0], np.array([(yy[i + 1] + yy[i]) / 2. for i in np.arange(ngauss_bins)])
def_xgauss = np.sum(
	np.exp(-xmid ** 2 / (2 * sigx ** 2)) * dx)  # for horizontal wires; amp = 1, apply real amplitude later
def_ygauss = np.sum(np.exp(-ymid ** 2 / (2 * sigy ** 2)) * dy)  # for vertical wires (unshadowed value)


def gauss_integral_check2():
	nchk = 777
	ww = np.linspace(-xgap / 2., xgap / 2., endpoint=True, num=nchk)
	dw, wmid = ww[1] - ww[0], [(ww[i + 1] + ww[i]) / 2. for i in np.arange(nchk - 1)]
	power2wire = [ygauss_amp(wmid[i]) * def_ygauss * dw for i in np.arange(nchk - 1)]  # integrated_power*area
	power = np.sum(power2wire)
	print(f'Ptot = {power * 1.e-3:.2f} (kW)')
	a = 1


# gauss_integral_check2()  # OK

vert_wire_power_computed = np.array([])
horz_wire_power_computed, power_shadowed = np.array([]), np.array([])
for wire_pos in xwire_pos_arr:  # compute power to vertical wires with NO shadowing
	ww = np.linspace(wire_pos - w_diam / 2., wire_pos + w_diam / 2., nwire_bins + 1)  # make nx bins btwn bounds of wire
	dw, wmid = ww[1] - ww[0], [(ww[i + 1] + ww[i]) / 2. for i in np.arange(nwire_bins)]
	power2wire = [ygauss_amp(wmid[i]) * def_ygauss * dw for i in np.arange(nwire_bins)]  # integrated_power*area
	vert_wire_power_computed = np.append(vert_wire_power_computed, np.sum(power2wire))
for wire_pos in ywire_pos_arr:  # compute power to horizontal wires
	ww = np.linspace(wire_pos - w_diam / 2., wire_pos + w_diam / 2., nwire_bins + 1)  # make nx bins btwn bounds of wire
	dw, wmid = ww[1] - ww[0], [(ww[i + 1] + ww[i]) / 2. for i in np.arange(nwire_bins)]
	power2wire = [xgauss_amp(wmid[i]) * def_xgauss * dw for i in np.arange(nwire_bins)]  # integrated_power*area
	horz_wire_power_computed = np.append(horz_wire_power_computed, np.sum(power2wire))
for wire_pos in xwire_pos_arr:  # compute shadowing per vertical wire
	pshadow = 0.  # shadowed power
	for ypos in ywire_pos_arr:
		pshadow += xygauss_amp(wire_pos, ypos) * w_diam ** 2  # assumes small wire diam so guass is ~const
	power_shadowed = np.append(power_shadowed, pshadow)
vert_wire_power_corrected = vert_wire_power_computed - power_shadowed

# compute predicted temperature change
# P/A = rho*l*cp*dT/tpulse from Titus' paper
# Q = m*c*dT
ywire_mass = np.pi * (w_diam / 2.) ** 2 * ygap * rho  # [g]
energy_to_ywires = vert_wire_power_corrected * tbeam  # [J]
dtempy = energy_to_ywires / ywire_mass / c_sp  # [K]
energy_to_ywires_uncorrected = vert_wire_power_computed * tbeam  # [J]
dtempy_uncorrected = energy_to_ywires_uncorrected / ywire_mass / c_sp  # [K]
xwire_mass = np.pi * (w_diam / 2.) ** 2 * xgap * rho
energy_to_xwires = horz_wire_power_computed * tbeam
dtempx = energy_to_xwires / xwire_mass / c_sp

# compute predicted max temp range (center of center wire)
# consider wire of length equal to its diameter (square area) at center of beam: power = A
energy_to_center = A * tbeam * w_diam * w_diam
mass_center = np.pi * (w_diam / 2.) ** 2 * w_diam * rho
dtemp_center = energy_to_center / mass_center / c_sp  # [K]

# compute predicted lengthening of wires due to increase temp
dlx = xgap * dtempx * alpha  # [cm]
dly = ygap * dtempy * alpha  # [cm]
dly_uncorr = ygap * dtempy_uncorrected * alpha  # [cm]

wire_intercept_fraction = np.sum([np.sum(horz_wire_power_computed), np.sum(vert_wire_power_corrected)]) / pbeam
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

if __name__ == '__main__':
	pass
