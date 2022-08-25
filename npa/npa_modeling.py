import datetime
from datetime import date

from transp_code.transp_classes import Halo3D, FBM
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker
from helpful_stuff import read_eqdsk, ltx_limiter, read_nenite
from npa.ionization_cross_sections import *
from scipy.interpolate import griddata, RegularGridInterpolator
import pickle
import os

'''
TODO:
define viewing aperture on detector view plot, then sum over window to get predicted flux and energy distribution
use NPA instrument function to predict current per detector
Use neutral density to calculate attenuation

THOUGHTS:
-f_reionized: Compute re-ionization fraction along path to detector (adjust # pts by distance to get a pt every cm or so)
 the same way I do in POET

DONE:
Include neutral density to compute emissivity
Include cross section data
Pull actual B field data (don't assume purely toroidal field)
Include viewing limitations: not all bmvol are able to emit onto detector
program in the sign convention for pitch
'''

# todo: use actual neutral halo data from TRANSP instead
# todo: NEED TO include rotation of transp halo data to account for generic placement of NPA
# WLOG set beam source at X=0, aiming defined by tangency radius
beamtan = 21.3  # [cm] (tangency radius of beam: RTCENA in TR.DAT)
beamsrc_to_tan = 257.  # [cm] (dist from beam source to tangency radius: XLBTNA in TR.DAT)
xsrc, ysrc = 0, np.sqrt(beamtan ** 2 + beamsrc_to_tan ** 2)
phi_src = np.arctan2(beamtan, beamsrc_to_tan)  # angle between src-machinecenter and beam-centerline


def get_neutral_density(x_bmvol, y_bmvol, z_bmvol, neut_halo):
	# mixing units here- all neut stuff in [cm] so convert bmvol from [m]->[cm]
	xbm, ybm, zbm = x_bmvol*1.e2, y_bmvol*1.e2, z_bmvol*1.e2
	r_pt2src = np.sqrt((xsrc - xbm) ** 2 + (ysrc - ybm) ** 2)
	phi_pt2src = np.arctan2(xsrc - xbm, ysrc - ybm)
	# compute transform of bmvol elements into box coords (x,y,z)->(x,l,y)
	# note bmvol_z = box_y
	# bb denotes .b.mvol element in .b.ox coords
	x_bb = r_pt2src * np.sin(phi_pt2src + phi_src)
	y_bb = zbm
	l_bb = r_pt2src * np.cos(phi_pt2src + phi_src)
	# return 1.e9  # [#/cm^3]
	if min(neut_halo.lbox) <= l_bb <= max(neut_halo.lbox) and min(neut_halo.ybox) <= y_bb <= max(
			neut_halo.ybox) and min(neut_halo.xbox) <= x_bb <= max(neut_halo.xbox):
		return fneut([l_bb, y_bb, x_bb])[0]
	else:
		return 1.  # sets background neutral density  (#/cm^3)


# Normal LTX ops is ip=cw, bt=ccw
ip_shouldbe_cw = True
bt_shouldbe_cw = False

# 100002F01 = 10ms 16keV 35A beam into 120kA Ip recon of 100981 @ 468ms eqdsk
# time points 1-5 at 0, 2.5, 5, 7.5, 10, 12.5ms (beam on at 0ms)
fbm_fn = '//samba/wcapecch/transp/t100002/100002F01_fi_4.cdf'
n0_fn = '//samba/wcapecch/transp/t103617/103617C01_boxn0_4.cdf'
neut_halo = Halo3D(n0_fn)
# only 1 beam box- remove extra dimension
neut_halo.boxn0 = neut_halo.boxn0[:, :, :, 0]
neut_halo.boxn0h0 = neut_halo.boxn0h0[:, :, :, 0]
neut_halo.boxn0hh = neut_halo.boxn0hh[:, :, :, 0]
neut_halo.total_neut = neut_halo.boxn0 + neut_halo.boxn0h0 + neut_halo.boxn0hh  # beam + 0 gen + higher gen neutrals (#/cm^3)
fneut = RegularGridInterpolator((neut_halo.lbox, neut_halo.ybox, neut_halo.xbox), neut_halo.total_neut)

eqdsk_fn = '//samba/wcapecch/datasets/LTX_100981_468-1_5.eqdsk'
eq = read_eqdsk(eqdsk_fn)
nenite_fn = '//samba/wcapecch/datasets/LTX_100981_468-1_5.nenite'
xb, plflx, ne_m3, ni_m3, te_ev = read_nenite(nenite_fn)
ne_m3, ni_m3 = ne_m3 * 1.e6, ni_m3 * 1.e6  # convert to m^-3
rlimiter, zlimiter, rminor_lim, theta_lim = ltx_limiter()
eq_r2d, eq_z2d, br2d, bz2d, bphi2d, psi_rz = eq['x_xy'], eq['y_xy'], eq['br_xy'], eq['bz_xy'], eq['bphi_xy'], eq[
	'psixy']
iz0 = np.where(abs(eq['y']) == min(abs(eq['y'])))[0][0]
psislice = psi_rz[iz0, :]  # take slice at midplane
if (np.where(psislice == min(psislice))[0][0] == 0) or (np.where(psislice == min(psislice))[0][0] == len(psislice)):
	psi_rz *= -1.  # if min is at endpoint, flip so min is in middle
psi_rz -= np.min(psi_rz)  # move offset to zero to match transp data

r_magax = 0.4  # doesn't need to be exact, just using interior pt. to find pts outside limiter
rminor_xy = np.sqrt((eq_r2d - r_magax) ** 2 + eq_z2d ** 2)
theta_xy = np.arctan2(eq_z2d, (eq_r2d - r_magax))

r_interp = np.interp(theta_xy, theta_lim, rminor_lim)  # rminor_limits of all theta pts of _xy
br2d, bz2d = np.ma.masked_where(rminor_xy > r_interp, br2d), np.ma.masked_where(rminor_xy > r_interp, bz2d)
bphi2d = np.ma.masked_where(rminor_xy > r_interp, bphi2d)

# field orientation handling
if eq['ip_is_cw'] != ip_shouldbe_cw:
	br2d *= -1.
	bz2d *= -1.
if eq['bt_is_cw'] != bt_shouldbe_cw:
	bphi2d *= -1.

rz_pts = np.append(eq_r2d.reshape((len(eq_r2d.flatten()), 1)), eq_z2d.reshape((len(eq_z2d.flatten()), 1)), axis=1)

magfig, (magax1, magax2, magax3) = plt.subplots(ncols=3, sharey=True)
magax1.contourf(eq_r2d, eq_z2d, np.ma.masked_where(rminor_xy > r_interp, br2d))
magax1.plot(rlimiter, zlimiter, 'k')
magax1.set_title('Br')
magax1.set_xlabel('R (m)')
magax1.set_ylabel('Z (m)')
magax2.contourf(eq_r2d, eq_z2d, np.ma.masked_where(rminor_xy > r_interp, bz2d))
magax2.plot(rlimiter, zlimiter, 'k')
magax2.set_title('Bz')
magax3.contourf(eq_r2d, eq_z2d, np.ma.masked_where(rminor_xy > r_interp, bphi2d))
magax3.plot(rlimiter, zlimiter, 'k')
magax3.set_title('B_tor')

# fbm has indices of [bmvol, pitch, energy]
if os.path.exists(fbm_fn.split('.')[0] + '.pkl'):  # pickled file exists, try loading this (might be faster)
	print('loading pickle file', end='')
	t1 = datetime.datetime.now()
	fbm = pickle.load(open(fbm_fn.split('.')[0] + '.pkl', 'rb'))
	t2 = datetime.datetime.now()
	print(f'... {(t2 - t1).seconds} sec')
else:
	print('no pickle file yet, creating one for later', end='')
	t1 = datetime.datetime.now()
	fbm = FBM(fbm_fn)
	pickle.dump(fbm, open(fbm_fn.split('.')[0] + '.pkl', 'wb'))
	t2 = datetime.datetime.now()
	print(f'... {(t2 - t1).seconds} sec')

# divide up toroidally (np.sum(bmvol) gives volume of torus)
nphi_tor = 100  # number of phi bins
phi_bin_limits = np.linspace(0, 2. * np.pi, num=nphi_tor + 1, endpoint=True)
dphi = 2. * np.pi / nphi_tor / 2.
phi_bins = phi_bin_limits[0:-1] + dphi
fbm.bmvol /= nphi_tor  # reduce bmvol volumes into nphi new segments

pitch_bin_limits = [(fbm.pitch[i] + fbm.pitch[i + 1]) / 2. for i in np.arange(len(fbm.pitch) - 1)]
pitch_bin_limits.insert(0, -1.)
pitch_bin_limits.append(1.)

mp = 1.672e-27  # proton mass [kg]
me = 9.109e-31  # electron mass [kg]
ee = 1.60218e-19  # charge [C]

emissiv = np.zeros_like(fbm.fbm[:, :, :])  # bmvol vs pitch, need energy axis since attenuation depends on energy
de_arr = [fbm.energy[e + 1] - fbm.energy[e] for e in np.arange(len(fbm.energy) - 1)]  # step size in energy array [eV]
if max(de_arr) - min(de_arr) > 1:  # diff in de > 1eV: stop
	raise ValueError('energy array on irregular grid')
de = np.mean(de_arr)
# cx and ion impact depend only on ion energy so we can pull these out and do them only once
sigv_cx_cm3_s = tabulated_i_ch_ex_ionization(fbm.energy, 1.)  # [cm^3/s]
# sigv_ii_cm3_s = tabulated_iimpact_ionization(fbm.energy, 1.)  # [cm^3/s]
print('computing emissivity of poloidal bmvols: num/sec when multiplied by n0 in #/cm^3')
for p in np.arange(len(fbm.pitch) - 1):
	for k in np.arange(len(fbm.bmvol)):
		# dtheta = abs(np.arccos(pitch_bin_limits[p]) - np.arccos(pitch_bin_limits[p + 1]))
		# theta_av = (np.arccos(pitch_bin_limits[p]) + np.arccos(pitch_bin_limits[p + 1])) / 2.
		# n0 = get_neutral_density(x_bmvol, y_bmvol, z_bmvol)  # [num/m^3] don't do this here- only doing cross section
		# v_elec = np.sqrt(2. * ee * te_fbmgrid[k] / me)
		# sigv_ei = tabulated_eimpact_ionization(fbm.energy, te_fbmgrid[k])  # [cm^3/s]
		v_ion_m_s = np.sqrt(2. * ee * fbm.energy / mp)  # m/s speed of ion
		# fbm has indices of [bmvol, pitch, energy]
		Ni_vs_energy = fbm.fbm[k, p, :] * fbm.bmvol[k] * de / len(fbm.pitch)  # [num at this pitch p in this bmvol k]
		num_sec_per_n0 = Ni_vs_energy * sigv_cx_cm3_s  # [num/sec when multiplied by n0 in #/cm^3]
		# a_annulus_type_thing = 2. * np.pi * dtheta * np.sin(theta_av)  # area of surf between pitch bin limits at R=1m
		# DON'T divide by annulus- keep as rate of emission per n0
		emissiv[k, p, :] = num_sec_per_n0  # / a_annulus_type_thing

# emissivity, det_phi, and det_theta for each bmvol element at each phi location
emissiv_det = np.zeros((len(fbm.bmvol), nphi_tor, len(fbm.energy)))
det_phi_theta_rmag = np.zeros((len(fbm.bmvol), nphi_tor, 3))

r_det, z_det = .7, 0.
x_det, y_det = 0., -r_det  # WLOG if we have toroidal symmetry
# get br, bz, bphi onto fbm grid
fbm_rz_pts = np.append(fbm.r2d.reshape((len(fbm.r2d), 1)), fbm.z2d.reshape((len(fbm.z2d), 1)), axis=1) * 1.e-2
br_fbmgrid = griddata(rz_pts, br2d.flatten(), fbm_rz_pts)
bz_fbmgrid = griddata(rz_pts, bz2d.flatten(), fbm_rz_pts)
bphi_fbmgrid = griddata(rz_pts, bphi2d.flatten(), fbm_rz_pts)
psi_fbmgrid = griddata(rz_pts, psi_rz.flatten(), fbm_rz_pts)
ne_fbmgrid_m3, ni_fbmgrid_m3, te_fbmgrid_ev = np.zeros_like(psi_fbmgrid), np.zeros_like(psi_fbmgrid), np.zeros_like(
	psi_fbmgrid)
for k in np.arange(len(psi_fbmgrid)):
	ne_fbmgrid_m3[k] = np.interp(psi_fbmgrid[k], plflx, ne_m3)
	ni_fbmgrid_m3[k] = np.interp(psi_fbmgrid[k], plflx, ni_m3)
	te_fbmgrid_ev[k] = np.interp(psi_fbmgrid[k], plflx, te_ev)

f = np.linspace(0, 1, endpoint=True)
print('computing emissivity [num/sec/m^2] onto detector per bin vs energy')
for iphi in np.arange(nphi_tor):
	print(f'{iphi + 1}'.zfill(2), end=' ' if (iphi + 1) % 25 != 0 else '\n')
	# print(f'{phi_bins[iphi]*180/np.pi:.5}', end=' ' if (iphi + 1) % 25 != 0 else '\n')
	for k in np.arange(len(fbm.bmvol)):
		r_bmvol, z_bmvol = fbm.r2d[k] * 1.e-2, fbm.z2d[k] * 1.e-2  # [m]
		x_bmvol, y_bmvol = -fbm.r2d[k] * 1.e-2 * np.sin(phi_bins[iphi]), fbm.r2d[k] * 1.e-2 * np.cos(
			phi_bins[iphi])  # [m]
		rx, ry, rz = x_det - x_bmvol, y_det - y_bmvol, z_det - z_bmvol
		r_mag = np.sqrt(rx ** 2 + ry ** 2 + rz ** 2)
		
		# check for visibility to detector
		r_line = np.sqrt((x_det + f * (x_bmvol - x_det)) ** 2 + (y_det + f * (y_bmvol - y_det)) ** 2)
		theta_line = np.arctan2(z_det + f * (z_bmvol - z_det), r_line - r_magax)
		rmin_line = np.sqrt((r_line - r_magax) ** 2 + (z_det + f * (z_bmvol - z_det)) ** 2)
		rminor_lim_interp = np.interp(theta_line, theta_lim, rminor_lim)
		if len(np.where(abs(np.diff(np.sign(rmin_line - rminor_lim_interp))) > 1)[0]) > 1:
			# multiple crossings of limiter (we allow one for detector)
			em_det = np.zeros_like(fbm.energy)
		else:
			# get B field
			bx = -br_fbmgrid[k] * np.sin(phi_bins[iphi]) - bphi_fbmgrid[k] * np.cos(phi_bins[iphi])
			by = br_fbmgrid[k] * np.cos(phi_bins[iphi]) - bphi_fbmgrid[k] * np.sin(phi_bins[iphi])
			bz = bz_fbmgrid[k]
			
			b_mag = np.sqrt(bx ** 2 + by ** 2 + bz ** 2)
			# TRANSP takes convention that pitch > 0 when aligned with PLASMA CURRENT
			if ip_shouldbe_cw != bt_shouldbe_cw:
				sign_convention = -1.  # since Bt is anti-aligned with LTX in general we apply a sign flip below
			else:
				sign_convention = 1.
			pitch_det = sign_convention * (rx * bx + ry * by + rz * bz) / r_mag / b_mag
			
			n0 = get_neutral_density(x_bmvol, y_bmvol, z_bmvol, neut_halo)  # [num/cm^3]
			emiss = emissiv[k, :, :] * n0  # [num/sec]
			em_det = np.zeros_like(fbm.energy)  # [num/sec/m^2]
			for ie in np.arange(len(fbm.energy)):
				em_det[ie] = np.interp(pitch_det, fbm.pitch,
				                       emiss[:, ie]) / r_mag ** 2  # flux falls of as 1/r^2 from bmvol element
		
		# detector viewing angles:
		det_phi = np.arctan(-rx / -ry)
		det_theta = np.arcsin(-rz / r_mag)
		emissiv_det[k, iphi, :] = em_det
		det_phi_theta_rmag[k, iphi, :] = [det_phi, det_theta, r_mag]

print('')
nph, nth = 50, 50
phi_det_bins = np.linspace(-np.pi / 2., np.pi / 2., endpoint=True, num=nph + 1)
phi_det_cntr = phi_det_bins[0:-1] + (phi_det_bins[1] - phi_det_bins[0]) / 2.
theta_det_bins = np.linspace(-np.pi, np.pi, endpoint=True, num=nth + 1)
theta_det_cntr = theta_det_bins[0:-1] + (theta_det_bins[1] - theta_det_bins[0]) / 2.
em2d = np.zeros((nph, nth))
print('analyzing bmvol elements from detector view')
for iph in np.arange(nph):  # emissiv_det = [bmvol, iphi, energy]
	for ith in np.arange(nth):
		# define path from det outward given angles phi, theta
		# get pts every cm or so along path from det out 2 meters
		r_ep = 2.  # set 2m length (long compared to LTX, ensuring we cover all bmvol elements)
		z_endpt = r_ep * np.sin(theta_det_cntr[ith])
		xy_endpt = r_ep * np.cos(theta_det_cntr[ith])
		x_endpt, y_endpt = xy_endpt * np.sin(phi_det_cntr[iph]), xy_endpt * np.cos(phi_det_cntr[iph])
		
		xpath = np.linspace(x_det, x_endpt, endpoint=False, num=int(r_ep * 100.))
		ypath = np.linspace(y_det, y_endpt, endpoint=False, num=int(r_ep * 100.))
		zpath = np.linspace(z_det, z_endpt, endpoint=False, num=int(r_ep * 100.))
		pts_along = np.linspace(0, r_ep, endpoint=False, num=int(r_ep * 100.))
		dx = np.sqrt((xpath[1] - xpath[0]) ** 2 + (ypath[1] - ypath[0]) ** 2 + (zpath[1] - zpath[0]) ** 2)
		rpath = np.sqrt(xpath ** 2 + ypath ** 2)
		phipath = np.arctan(-xpath, ypath)
		rz_path = np.append(rpath.reshape((len(rpath), 1)), zpath.reshape((len(zpath), 1)), axis=1)
		psi_path = griddata(rz_pts, psi_rz.flatten(), rz_path)
		ni_path_m3 = np.zeros_like(psi_path)
		# ne_path_m3, te_path_ev = np.zeros_like(psi_path), np.zeros_like(psi_path)
		for ipt in np.arange(len(psi_path)):
			# ne_path_m3[ipt] = np.interp(psi_path[ipt], plflx, ne_m3)
			ni_path_m3[ipt] = np.interp(psi_path[ipt], plflx, ni_m3)
		# te_path_ev[ipt] = np.interp(psi_path[ipt], plflx, te_ev)
		# sigv_ei = tabulated_eimpact_ionization(fbm.energy, te_path)
		frac_ioniz = np.outer(ni_path_m3, (sigv_cx_cm3_s * 1.e-6) / v_ion_m_s * dx)  # vs position and energy
		frac_ioniz = np.nan_to_num(frac_ioniz)  # sets nan to 0.
		frac_trans = 1. - frac_ioniz
		frac_trans = np.cumprod(frac_trans, axis=0)  # cumulative product along position vector v. energy
		
		iwin = np.where(
			(det_phi_theta_rmag[:, :, 0] >= phi_det_bins[iph]) & \
			(det_phi_theta_rmag[:, :, 0] < phi_det_bins[iph + 1]) & \
			(det_phi_theta_rmag[:, :, 1] >= theta_det_bins[ith]) & \
			(det_phi_theta_rmag[:, :, 1] < theta_det_bins[ith + 1]))
		
		em_sum = 0.
		for ibmvol in np.arange(len(iwin[0])):
			emiss_v_energy = emissiv_det[iwin[0][ibmvol], iwin[1][ibmvol], :]
			r_to_bmvol = det_phi_theta_rmag[iwin[0][ibmvol], iwin[1][ibmvol], 2]
			iclos = np.where(abs(pts_along - r_to_bmvol) == min(abs(pts_along - r_to_bmvol)))[0][0]
			em_sum += np.sum(frac_trans[iclos, :] * emiss_v_energy)
		em2d[iph, ith] = em_sum

print('plotting')
# detector view
fig, ax = plt.subplots()
em2d = np.ma.masked_where(em2d <= 0, em2d)
cb1 = ax.contourf(phi_det_cntr, theta_det_cntr, em2d.transpose(), locator=ticker.LogLocator())
cbar1 = fig.colorbar(cb1)
cbar1.set_label('num/sec/m^2')
ax.axvline(0, c='k', ls='--', linewidth=1)
ax.axhline(0, c='k', ls='--', linewidth=1)
# ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
# ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.set_xlabel('$\\phi$')
ax.set_ylabel('$\\theta$')

nr = 50
rr_bins = np.linspace(0, .7, endpoint=True, num=nr + 1)
rr = np.linspace(0, .7, endpoint=False, num=nr) + .7 / 50 / 2.
r2, p2 = np.meshgrid(rr, phi_bin_limits + dphi)

# top down view of emissivity summed over all pitch
# fig3, ax3 = plt.subplots(subplot_kw=dict(projection='polar'))
# e1 = np.zeros((nphi_tor, len(rr)))
# for iph in np.arange(nphi_tor):
# 	for ir in np.arange(len(rr)):
# 		e1[iph, ir] = np.sum(emissiv[np.where(
# 			(fbm.r2d * 1.e-2 > rr_bins[ir]) & (fbm.r2d * 1.e-2 < rr_bins[ir + 1])
# 		), :])
# e1 = np.vstack([e1, e1[0, :]])
# e1 = np.ma.masked_where(e1 <= 0, e1)
# cb3 = ax3.contourf(p2, r2, e1, locator=ticker.LogLocator())
# cbar3 = fig3.colorbar(cb3)
# cbar3.set_label('num/sec')

# plot of 3d Halo Neutral Density
fig3, ax3 = plt.subplots(subplot_kw=dict(projection='polar'))
n0tot = neut_halo.total_neut[:, :, :]
n0tot_mp = np.sum(n0tot, axis=1)  # sum over vertical coord
xbox, lbox = neut_halo.xbox, neut_halo.lbox
x_bb, y_bb = np.meshgrid(xbox, lbox)
r_bb, phi_bb = np.zeros_like(x_bb), np.zeros_like(x_bb)
for ix in np.arange(len(xbox)):
	for il in np.arange(len(lbox)):
		x_bb[il, ix] = xsrc + lbox[il] * np.sin(phi_src) - xbox[ix] * np.cos(phi_src)
		y_bb[il, ix] = ysrc - lbox[il] * np.cos(phi_src) - xbox[ix] * np.sin(phi_src)
		r_bb[il, ix] = np.sqrt(x_bb[il, ix] ** 2 + y_bb[il, ix] ** 2)
		phi_bb[il, ix] = np.arctan2(-x_bb[il, ix], y_bb[il, ix])  # neg here to get NB injecting co-Ip
# all neut stuff in [cm], convert to [m] to be consistent w/other plot
r_bb *= 1.e-2
cb3 = ax3.contourf(phi_bb, r_bb, n0tot_mp, levels=25)
cbar3 = fig3.colorbar(cb3)
cbar3.set_label('#/cm^3 (summed over vertical coord)')
ax3.set_rlim((0, .75))
ax3.set_theta_zero_location('N')
ax3.set_title('beam neutral halo')

# top down view of emissivity reaching detector
em_top = np.zeros((nphi_tor + 1, len(rr)))
for iph in np.arange(nphi_tor):
	for ir in np.arange(len(rr)):
		em_top[iph, ir] = np.sum(emissiv_det[:, iph, 0][np.where(
			(fbm.r2d * 1.e-2 > rr_bins[ir]) & (fbm.r2d * 1.e-2 < rr_bins[ir + 1])
		)])
		if iph == 0:
			em_top[nphi_tor, ir] = em_top[iph, ir]
# em_top = np.vstack([em_top, em_top[0, :]])
fig2, ax2 = plt.subplots(subplot_kw=dict(projection='polar'))
phi_det = np.arctan2(-x_det, y_det)
em_top = np.ma.masked_where(em_top <= 0, em_top)
cb2 = ax2.contourf(p2, r2, em_top, locator=ticker.LogLocator())
cbar2 = fig2.colorbar(cb2)
cbar2.set_label('num/sec/m^2')
ax2.plot(phi_det, r_det, 'ro')
ax2.annotate('detector', (phi_det, r_det), xytext=(.7, -.9), textcoords='figure fraction',
             arrowprops=dict(facecolor='r', shrink=.05))
ax2.set_theta_zero_location('N')

plt.show()
