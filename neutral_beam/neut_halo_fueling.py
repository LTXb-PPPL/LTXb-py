import pickle
import random

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from toolbox.ionization_cross_sections import *
from toolbox.helpful_stuff import read_nenite, ltx_limiter, read_eqdsk3, closest, ltx_above
from toolbox.particle_trajectory_tools import get_generic_path
from transp_code.transp_classes import Halo3D

"""
We have neutral footprint from beam via Halo model in NUBEAM
TRANSP monitors neutral tracks, and so must keep track of ionization events, but does not track fueling vs non-fueling
type ionizations.
So the idea here is to take the neutral halo and compute the following fractions:
flost = fraction lost via collision with boundary
ffuel = fraction that eventually ionized via electron impact (ion impact is negligible at thermal temps)
When the beam enters the torus, a certain fraction ionizes either through CX or IE/II events. The IE/II fraction
(roughly 1/3) contributes directly to fueling, and the ffuel fraction determined here computes what fraction of the CX
ionizations ultimately lead to a fueling event.
So ultimately the beam fueling fraction will be:
beam_fueling_fraction = frac_initial_cx*ffuel + frac_initial_impact_ionization
"""

'''CONTROLS'''
new2dscan = True
newhalo = True
nth = 1  # downsample to every nth point
tol = 1.e-4  # tolerance for convergence of neutral tracks at given location
''' VARIABLES '''
shot = 106536
direc = 'Z:/transp/'
# fsav2d = f'Z:/PycharmProjects/LTXb-py/neutral_beam/neut_halo_fueling_ds{nth}.pkl'
fsav2d = f'Z:/PycharmProjects/LTXb-py/neutral_beam/neut_halo_fueling_ds{nth}_10xdensity.pkl'
for (rtanfn, rtan) in zip([2], [24.]):
	fbm_fn = f'{direc}t{shot}/{shot}R02_fi_5.cdf'
	n0_fn = f'{direc}t{shot}/{shot}R02_boxn0_5.cdf'
	eq_fn = f'{direc}t{shot}/{shot}R02_05.eqdsk'
	nenite_fn = f'{direc}t{shot}/{shot}R02_05.nenite'
ti = 100.  # eV temp of thermal ions
te = 150.  # eV temp of thermal electrons
ion_amu = 1.  # amu of ion species
'''***************'''

clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']

# interp density onto equilibrium poloidal plane
eq = read_eqdsk3(eq_fn)
xb, plflx, ne_cm3, ni_cm3, te_ev = read_nenite(nenite_fn)
# artificially play w/density
ne_cm3 *= 10.
ni_cm3 *= 10.
# ne_m3, ni_m3 = ne_cm3 * 1.e6, ni_cm3 * 1.e6  # convert to m^-3
rlimiter, zlimiter, rminor_lim, theta_lim = ltx_limiter()
eq_r2d, eq_z2d, br2d, bz2d, bphi2d, psi_rz = eq['x_xy'], eq['y_xy'], eq['br_xy'], eq['bz_xy'], eq['bphi_xy'], eq[
	'psixy']
r1d, z1d = eq['x'], eq['y']

# reduce grid resolution for faster model
eq_r2d, eq_z2d, br2d = eq_r2d[::nth, ::nth], eq_z2d[::nth, ::nth], br2d[::nth, ::nth]
bz2d, bphi2d, psi_rz = bz2d[::nth, ::nth], bphi2d[::nth, ::nth], psi_rz[::nth, ::nth]
r1d, z1d = r1d[::nth], z1d[::nth]

ne_2d, ni_2d, fuel_2d = np.zeros_like(eq_r2d), np.zeros_like(eq_r2d), np.zeros_like(eq_r2d)
for ir in np.arange(len(r1d)):
	for iz in np.arange(len(z1d)):
		ne_2d[iz, ir] = np.interp(psi_rz[iz, ir], plflx, ne_cm3)
		ni_2d[iz, ir] = np.interp(psi_rz[iz, ir], plflx, ni_cm3)
ni_func = RegularGridInterpolator((z1d, r1d), ni_2d)  # functions to pull n values from later
ne_func = RegularGridInterpolator((z1d, r1d), ne_2d)

# compute mfp over poloidal cross section
sigcx = tabulated_i_ch_ex_ionization(ti, ion_amu, retrn='sigma')  # cm^2
sigvei = tabulated_eimpact_ionization(ti, te)  # cm^3/s
vel_i = np.sqrt(2. * qe * ti / (mi * ion_amu)) * 100.  # cm/s
sigei = sigvei / vel_i

cx_mfp = 1. / (sigcx * ni_2d)
ei_mfp = 1. / (sigei * ne_2d)  # electron-ion impact mean free path
sigtot = sigcx + sigei
mfptot = 1. / (sigcx * ni_2d + sigei * ne_2d)
ei_frac = mfptot / ei_mfp  # (num/m ei events) / (num/m total events)


def pt_inside_boundary(r, z):
	rmag = 0.4  # doesn't need to be exact, just using interior pt. to find pts outside limiter
	rmin = np.sqrt((r - rmag) ** 2 + z ** 2)
	thet = np.arctan2(z, (r - rmag))
	if rmin < np.interp(thet, theta_lim, rminor_lim):
		inside = True
	else:
		inside = False
	return inside


def truncate_at_boundary(rr, zz):
	rmag = 0.4
	rmin = np.sqrt((rr - rmag) ** 2 + zz ** 2)
	thet = np.arctan2(zz, (rr - rmag))
	rminor_lim_interp = np.interp(thet, theta_lim, rminor_lim)  # interpolated values of limiter rmin on thet array
	itrunc = np.where(rmin > rminor_lim_interp)[0]
	if len(itrunc) > 0:  # some points outside boundary- should always be the case
		return itrunc[0]  # return index of first pt outside boundary
	else:
		raise ValueError('Path was not truncated- it should be truncated. Suggest investigating.')


def compute_neut_track(r, z):
	particle_fraction_fueled = 0.
	particle_fraction, min_particle_frac = 1., 1.e-3
	while particle_fraction > min_particle_frac:  # loop until the fraction of particle we're tracking is very small
		# print(f'particle_fraction: {particle_fraction}')
		xp, yp, rp, zp, dp = get_generic_path(r, 0., z)  # WLOG set x=r, y=0 (axisymmetry); dp=path step size
		itrunc = truncate_at_boundary(rp, zp)
		if itrunc == 0:
			particle_fraction = 0.  # path outside
			break
		xp, yp, rp, zp = xp[:itrunc - 1], yp[:itrunc - 1], rp[:itrunc - 1], zp[:itrunc - 1]
		visualize_path = False
		if visualize_path:
			xi, xo, yi, yo = ltx_above()
			fig, (ax1, ax2) = plt.subplots(ncols=2)
			ax1.plot(rlimiter, zlimiter)
			ax1.plot(rp, zp)
			ax2.plot(xi, yi)
			ax2.plot(xo, yo)
			ax2.plot(xp, yp)
			ax1.plot(rp[:itrunc - 1], zp[:itrunc - 1], 'r-')
		good_density = False
		while not good_density:  # don't know why but these density calls fail occasionally but just need to be rerun?
			try:
				n_i = np.array([ni_func([zp[i], rp[i]]) for i in range(len(rp))])
				n_e = np.array([ne_func([zp[i], rp[i]]) for i in range(len(rp))])
				good_density = True
			except ValueError:
				good_density = False
		local_frac_ionized = (n_i * sigcx + n_e * sigei) * dp
		frac_lost = 1. - np.sum(local_frac_ionized)  # ignore fraction lost
		frac_cx = np.sum(n_i * sigcx * dp)  # this becomes portion of particle to continue modeling
		frac_ei = np.sum(n_e * sigei * dp)  # this is portion of particle that fuels
		particle_fraction_fueled += frac_ei * particle_fraction
		particle_fraction = frac_cx * particle_fraction  # now only follow portion of particle that CX
	return particle_fraction_fueled


visualize_grid = 0
if visualize_grid:
	plt.plot(eq_r2d, eq_z2d, 'r.')
	plt.plot(rlimiter, zlimiter, 'k-')

# Go through each R, Z pt and compute loss/fuel fractions
if new2dscan:
	for ir in range(len(r1d)):
		for iz in range(len(z1d)):
			# ir, iz = 4, 39  # FOR TESTING
			print(f'{ir}/{len(r1d)}:{iz}/{len(z1d)}')
			ok = pt_inside_boundary(r1d[ir], z1d[iz])
			if z1d[iz] < 0:  # consider only upper half of poloidal cross section (up-down symmetry)
				ok = False
			if ok:
				ffuel_rz = [0.]
				eps = 1.
				ncomp, nfuel = 0, 0
				while eps > tol or ncomp < 20:
					for ii in range(5):  # compute 5 more tracks
						neut_fueled = compute_neut_track(r1d[ir], z1d[iz])
						ncomp += 1
						nfuel += neut_fueled  # will be fraction- following portion of particles that fuel
					# print(f'{ncomp} tracks computed: eps={eps}')
					eps = abs(ffuel_rz[-1] - nfuel / ncomp)  # how much has fueling fraction changed
					ffuel_rz.append(nfuel / ncomp)  # update fueling fraction
				fuel_2d[iz, ir] = ffuel_rz[-1]
			else:  # started outside boundary
				fuel_2d[iz, ir] = np.nan
	
	# since we only looked at area above midplane, mirror below for full poloidal cross section
	for i in range(len(z1d)):
		if z1d[i] < 0:
			fuel_2d[i, :] = fuel_2d[-1 - i, :]
	
	dat2d = {'r2d': eq_r2d, 'z2d': eq_z2d, 'fuel2d': fuel_2d}
	pickle.dump(dat2d, open(fsav2d, 'wb'))
	print(f'saved: {fsav2d}')
else:
	dat2d = pickle.load(open(fsav2d, 'rb'))
	eq_r2d, eq_z2d, fuel_2d = dat2d['r2d'], dat2d['z2d'], dat2d['fuel2d']

if newhalo:
	# SIMILAR to npa_modeling.py handling, create fneut function to compute neut density as function of
	halo = Halo3D(n0_fn)  # [NLBOX, NYBOX, NXBOX, NBBOX]
	# only 1 beam box- remove extra dimension
	halo.boxn0 = halo.boxn0[:, :, :, 0]
	halo.boxn0h0 = halo.boxn0h0[:, :, :, 0]
	halo.boxn0hh = halo.boxn0hh[:, :, :, 0]
	halo.total_neut = halo.boxn0 + halo.boxn0h0 + halo.boxn0hh  # beam + 0 gen + higher gen neutrals (#/cm^3)
	halo.thermal_neut = halo.boxn0h0 + halo.boxn0hh  # ignore beam neutrals (high energy), only look at secondaries
	
	beamtan = 24.  # [cm] RTCENA(1) in TR.DAT
	beamsrc_to_tan = 257.  # [cm] (dist from beam source to tangency radius: XLBTNA in TR.DAT)
	xsrc, ysrc = 0, np.sqrt(beamtan ** 2 + beamsrc_to_tan ** 2)
	phi_src = np.arctan2(beamtan, beamsrc_to_tan)  # angle between src-machinecenter and beam-centerline
	
	halo2d = np.zeros_like(ne_2d)
	for il in range(len(halo.lbox)):
		for ix in range(len(halo.xbox)):
			for iy in range(len(halo.ybox)):
				xpt = xsrc + halo.lbox[il] * np.sin(phi_src) - halo.xbox[ix] * np.cos(phi_src)
				ypt = ysrc - halo.lbox[il] * np.cos(phi_src) - halo.xbox[ix] * np.sin(phi_src)
				rpt = np.sqrt(xpt ** 2 + ypt ** 2) * 1.e-2
				zpt = halo.ybox[iy] * 1.e-2  # convert to m to be consistent w/r1d,z1d
				npt = halo.thermal_neut[il, iy, ix]  # num/cm^3
				if pt_inside_boundary(rpt, zpt):
					(_, ir) = closest(r1d, rpt)
					(_, iz) = closest(z1d, zpt)
					halo2d[iz, ir] += npt  # add neut density to 2d slice
	dv = halo.dv
	dat2d = {'r2d': eq_r2d, 'z2d': eq_z2d, 'fuel2d': fuel_2d, 'halo2d': halo2d, 'dv': dv}
	pickle.dump(dat2d, open(fsav2d, 'wb'))
	print(f'saved: {fsav2d}')
else:
	dat2d = pickle.load(open(fsav2d, 'rb'))
	eq_r2d, eq_z2d, fuel_2d, halo2d, dv = dat2d['r2d'], dat2d['z2d'], dat2d['fuel2d'], dat2d['halo2d'], dat2d['dv']

ntot = np.sum(halo2d.flatten()) * dv  # total neuts in halo
neut_fueled = fuel_2d / 100. * halo2d  # num/cm^3
ntot_fuel = np.nansum(neut_fueled.flatten()) * dv  # total neuts fueled
print(f'Ntot: {ntot}, Nfuel: {ntot_fuel} ({ntot_fuel / ntot * 100.}%)')

fig, ((ax2, ax3), (ax1, ax4)) = plt.subplots(nrows=2, ncols=2, sharey='row', sharex='col')
ax4.axis('off')
c = ax1.contourf(eq_r2d, eq_z2d, fuel_2d * 100.)
cb = fig.colorbar(c, ax=ax1)
cb.set_label('Fueling (%)')
# fig, ax2 = plt.subplots()
c2 = ax2.contourf(eq_r2d, eq_z2d, halo2d)
cb2 = fig.colorbar(c2, ax=ax2)
cb2.set_label('Halo Neutral Density ($cm^-3$)')
c3 = ax3.contourf(eq_r2d, eq_z2d, neut_fueled)
cb3 = fig.colorbar(c3, ax=ax3)
cb3.set_label('Halo fueling ($cm^-3$)')
for ax in [ax1, ax2, ax3]:
	# ax.axis('equal')
	ax.set_xlabel('R (m)')
	ax.set_ylabel('Z (m)')
	ax.set_ylim((-.4, .4))

plt.tight_layout()
plt.show()

if __name__ == '__main__':
	pass
