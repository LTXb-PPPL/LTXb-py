import matplotlib.pyplot as plt
import numpy as np

'''
given beam footprint (assumed gaussian), compute power intercepted by wire at various positions
generalize to independent x, y gaussian
This is to model predicted dT for given wire dimensions
'''
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']

# controls
nwire_bins = 100
ngauss_bins = 100
do_convergence_check = False

# beam characteristics
fwhmx = 8.  # cm
fwhmy = 9.  # cm
sigx = fwhmx / (2. * np.sqrt(2. * np.log(2)))  # [cm]
sigy = fwhmy / (2. * np.sqrt(2. * np.log(2)))  # [cm]
pbeam = 400.e3  # [W] total beam power
tbeam = 10.e-3  # [s] beam pulse duration

# wire characteristics
w_wid, w_len = .05, 32  # [cm] TAE design, 5mm diam wires 32 cm long
xw = -2.5  # x loc of wire

# 2d gaussian form- no covariance
# z = A * np.exp(-(x - xw) ** 2 / (2. * sigx ** 2) - (y - yw) ** 2 / (2. * sigy ** 2))
# so for given values of x, xw, sigx then z becomes 1d gaussian in y
# total integral: V=2*pi*A*sigx*sigy=total beam power
A = pbeam / (2. * np.pi * sigx * sigy)


def ygauss_amp(x):
	# return amplitude of guassian in y given location x
	return A * np.exp(-(x) ** 2 / (2. * sigx ** 2))


def ygauss_int(x):
	# return indefinite integral of gaussian in y at location x
	return A * np.exp(-x ** 2 / (2. * sigx ** 2)) * np.sqrt(2. * np.pi * sigy ** 2)


if do_convergence_check:
	fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
	# number of bins across the wire
	nxarr = np.arange(10, 210, 10)
	amp = ygauss_amp(xw)
	power_computed = []
	gauss_integral = []
	for nx in nxarr:
		xx = np.linspace(xw - w_wid / 2., xw + w_wid / 2., nx + 1)  # make nx bins between bounds of wire
		dx, xmid = xx[1] - xx[0], [(xx[i + 1] + xx[i]) / 2. for i in np.arange(nx)]
		power2wire = [ygauss_int(xmid[i]) * dx for i in np.arange(nx)]
		power_computed.append(np.sum(power2wire))
		yy = np.linspace(-w_len / 2., w_len / 2., endpoint=True, num=nx + 1)
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
yy = np.linspace(-w_len / 2., w_len / 2., endpoint=True, num=ngauss_bins + 1)
dy, ymid = yy[1] - yy[0], np.array([(yy[i + 1] + yy[i]) / 2. for i in np.arange(ngauss_bins)])
def_ygauss = np.sum(
	np.exp(-ymid ** 2 / (2 * sigy ** 2)))  # for vertical wires; amplitude = 1 here, apply real amplitude later
def_xgauss = np.sum(np.exp(-ymid ** 2 / (2 * sigx ** 2)))  # for horizontal wires

wire_pos_arr = np.arange(-15, 20, 2.5)  # every 5cm out to 25cm
power_computed = np.array([])
for wire_pos in wire_pos_arr:
	xx = np.linspace(wire_pos - w_wid / 2., wire_pos + w_wid / 2., nwire_bins + 1)  # make nx bins between bounds of wire
	dx, xmid = xx[1] - xx[0], [(xx[i + 1] + xx[i]) / 2. for i in np.arange(nwire_bins)]
	power2wire = [ygauss_amp(xmid[i]) * def_ygauss * dx for i in np.arange(nwire_bins)]
	power_computed = np.append(power_computed, np.sum(power2wire))

# compute predicted temperature change
# P/A = rho*l*cp*dt/tpulse
w_area = w_len*w_wid*1.e-4  # [m^2]
dtemp = pbeam*tbeam/w_area/

plt.plot(wire_pos_arr, power_computed * 1.e-3, 'o-')
plt.xlabel('wire position (cm)')
plt.ylabel('power to wire (kW)')
plt.show()

if __name__ == '__main__':
	pass
