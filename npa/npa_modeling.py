import datetime
from datetime import date

from npa.npa_channel_efficiencies import NPA
from transp_code.transp_classes import Halo3D, FBM
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker
from helpful_stuff import read_eqdsk3, ltx_limiter, read_nenite
from npa.ionization_cross_sections import *
from scipy.interpolate import griddata, RegularGridInterpolator
import pickle
import os

# todo: why is em2d empty with phi_det_degree = 170.?
'''
TODO:
Compute fast ion distribution and neutral halo for rtan=21.3, 24, 33cm to compare energy distribution of NPA signal
ALL emissivity is evaluated over grid in theta/phi or x/y- need to convert num/sec/m^2 to num/sec/m^2/sr
	-should be easy on detector view plot: dtheta and dphi map to a specific sr value
	-similar on top view plot, we're scanning in phi/theta
I think this is only CX- incorporate II and IE impact ionization
define viewing aperture on detector view plot, then sum over window to get predicted flux and energy distribution
use NPA instrument function to predict current per detector

THOUGHTS:

DONE:
f_reionized: Compute re-ionization fraction along path to detector
Use neutral density to calculate attenuation (?)
Include neutral density to compute emissivity
Include cross section data
Pull actual B field data (don't assume purely toroidal field)
Include viewing limitations: not all bmvol are able to emit onto detector
program in the sign convention for pitch
'''


def synthetic_npa(fbm_fn, n0_fn, eq_fn, nenite_fn, aperture_view=None, savfn='Z:/default.pkl', beamtan=21.3):
	# NPA detector position
	tangential = 0
	if tangential:
		r_det, phi_det_degree, z_det = .8, 190., 0.  # tangential view
		# aperture_phi0_deg, aperture_theta0_deg = 15, 0  # rtan ~ 20
		aperture_phi0_deg, aperture_theta0_deg = 24, 0  # rtan ~ 33
	else:
		r_det, phi_det_degree, z_det = .4, 225., .8  # tangential view
		aperture_phi0_deg, aperture_theta0_deg = 0, -90  # view down
	
	aperture_phi0, aperture_theta0 = aperture_phi0_deg * np.pi / 180., aperture_theta0_deg * np.pi / 180.
	aperture_view = [aperture_phi0, aperture_theta0]
	
	# r_det, phi_det_degree, z_det = .5, 190., .4
	x_det, y_det = -r_det * np.sin(phi_det_degree * np.pi / 180.), r_det * np.cos(phi_det_degree * np.pi / 180.)
	phi_det2machinecenter = np.arctan2(-x_det, y_det) + np.pi  # angle from det to machine center
	
	# Normal LTX ops is ip=cw, bt=ccw
	ip_shouldbe_cw = True
	bt_shouldbe_cw = False
	# WLOG set beam source at X=0, aiming defined by tangency radius
	# default: beamtan = 21.3  # [cm] (tangency radius of beam: RTCENA in TR.DAT)
	beamsrc_to_tan = 257.  # [cm] (dist from beam source to tangency radius: XLBTNA in TR.DAT)
	xsrc, ysrc = 0, np.sqrt(beamtan ** 2 + beamsrc_to_tan ** 2)
	phi_src = np.arctan2(beamtan, beamsrc_to_tan)  # angle between src-machinecenter and beam-centerline
	
	if not savfn.endswith('.pkl'):
		savfn = f'{savfn}.pkl'
	
	if not os.path.exists(savfn):
		neut_halo = Halo3D(n0_fn)
		# only 1 beam box- remove extra dimension
		neut_halo.boxn0 = neut_halo.boxn0[:, :, :, 0]
		neut_halo.boxn0h0 = neut_halo.boxn0h0[:, :, :, 0]
		neut_halo.boxn0hh = neut_halo.boxn0hh[:, :, :, 0]
		neut_halo.total_neut = neut_halo.boxn0 + neut_halo.boxn0h0 + neut_halo.boxn0hh  # beam + 0 gen + higher gen neutrals (#/cm^3)
		fneut = RegularGridInterpolator((neut_halo.lbox, neut_halo.ybox, neut_halo.xbox), neut_halo.total_neut)
		
		def get_neutral_density(x_bmvol, y_bmvol, z_bmvol, neut_halo):
			# mixing units here- all neut stuff in [cm] so convert bmvol from [m]->[cm]
			xbm, ybm, zbm = x_bmvol * 1.e2, y_bmvol * 1.e2, z_bmvol * 1.e2
			r_pt2src = np.sqrt((xsrc - xbm) ** 2 + (ysrc - ybm) ** 2)
			phi_pt2src = np.arctan2(xsrc - xbm, ysrc - ybm)
			# compute transform of bmvol elements into box coords (x,y,z)->(x,l,y)
			# note bmvol_z = box_y
			# bb denotes .b.mvol element in .b.ox coords
			x_bb = r_pt2src * np.sin(phi_pt2src + phi_src)
			y_bb = zbm
			l_bb = r_pt2src * np.cos(phi_pt2src + phi_src)
			if min(neut_halo.lbox) <= l_bb <= max(neut_halo.lbox) and min(neut_halo.ybox) <= y_bb <= max(
					neut_halo.ybox) and min(neut_halo.xbox) <= x_bb <= max(neut_halo.xbox):
				return fneut([l_bb, y_bb, x_bb])[0]
			else:
				return 1.e9  # sets background neutral density  (#/cm^3)
		
		eq = read_eqdsk3(eq_fn)
		xb, plflx, ne_m3, ni_m3, te_ev = read_nenite(nenite_fn)
		ne_m3, ni_m3 = ne_m3 * 1.e6, ni_m3 * 1.e6  # convert to m^-3
		rlimiter, zlimiter, rminor_lim, theta_lim = ltx_limiter()
		eq_r2d, eq_z2d, br2d, bz2d, bphi2d, psi_rz = eq['x_xy'], eq['y_xy'], eq['br_xy'], eq['bz_xy'], eq['bphi_xy'], \
		                                             eq[
			                                             'psixy']
		iz0 = np.where(abs(eq['y']) == min(abs(eq['y'])))[0][0]
		psislice = psi_rz[iz0, :]  # take slice at midplane
		if (np.where(psislice == min(psislice))[0][0] == 0) or (
				np.where(psislice == min(psislice))[0][0] == len(psislice)):
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
		
		rz_pts = np.append(eq_r2d.reshape((len(eq_r2d.flatten()), 1)), eq_z2d.reshape((len(eq_z2d.flatten()), 1)),
		                   axis=1)
		
		plot_mag = False
		if plot_mag:
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
		fbm = FBM(fbm_fn)
		
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
		
		emissiv = np.zeros_like(
			fbm.fbm[:, :, :])  # bmvol vs pitch, need energy axis since attenuation depends on energy
		de_arr = [fbm.energy[e + 1] - fbm.energy[e] for e in
		          np.arange(len(fbm.energy) - 1)]  # step size in energy array [eV]
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
				Ni_vs_energy = fbm.fbm[k, p, :] * fbm.bmvol[k] * de / len(
					fbm.pitch)  # [num at this pitch p in this bmvol k]
				num_sec_per_n0 = Ni_vs_energy * sigv_cx_cm3_s  # [num/sec when multiplied by n0 in #/cm^3]
				emissiv[k, p, :] = num_sec_per_n0
		
		emissiv_det = np.zeros((len(fbm.bmvol), nphi_tor, len(fbm.energy)))
		det_phi_theta_rmag = np.zeros((len(fbm.bmvol), nphi_tor, 3))
		
		# get br, bz, bphi onto fbm grid
		fbm_rz_pts = np.append(fbm.r2d.reshape((len(fbm.r2d), 1)), fbm.z2d.reshape((len(fbm.z2d), 1)), axis=1) * 1.e-2
		br_fbmgrid = griddata(rz_pts, br2d.flatten(), fbm_rz_pts)
		bz_fbmgrid = griddata(rz_pts, bz2d.flatten(), fbm_rz_pts)
		bphi_fbmgrid = griddata(rz_pts, bphi2d.flatten(), fbm_rz_pts)
		psi_fbmgrid = griddata(rz_pts, psi_rz.flatten(), fbm_rz_pts)
		ne_fbmgrid_m3, ni_fbmgrid_m3, te_fbmgrid_ev = np.zeros_like(psi_fbmgrid), np.zeros_like(
			psi_fbmgrid), np.zeros_like(
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
				
				# viewing angles from detector to bmvol element:
				det_phi = np.arctan2(-rx, -ry)  # angle from det to bmvol element using phi=0=North coords
				det_phi += phi_det2machinecenter  # angle from det to bmvol element relative to machinecenter
				det_theta = np.arcsin(-rz / r_mag)  # VERTICAL ANGLE
				emissiv_det[k, iphi, :] = em_det
				det_phi_theta_rmag[k, iphi, :] = [det_phi, det_theta, r_mag]
		
		# todo: determine min/max of phi and theta ranges based on bmvol angles to detector- emcompass all but no extra space
		print('analyzing all bmvol elements from detector view', end='')
		t1 = datetime.datetime.now()
		nph, nth = 51, 25
		phi_det_bins = np.linspace(-np.pi, np.pi, endpoint=True, num=nph + 1)  # MIDPLANE ANGLE
		phi_det_cntr = phi_det_bins[0:-1] + (phi_det_bins[1] - phi_det_bins[0]) / 2.
		theta_det_bins = np.linspace(-np.pi / 2., np.pi / 2., endpoint=True, num=nth + 1)  # VERTICAL ANGLE
		theta_det_cntr = theta_det_bins[0:-1] + (theta_det_bins[1] - theta_det_bins[0]) / 2.
		em2d = np.zeros((nph, nth))
		for iph in np.arange(nph):  # emissiv_det = [bmvol, iphi, energy]
			for ith in np.arange(nth):
				# define path from det outward given angles phi, theta
				# get pts every cm or so along path from det out 2 meters
				r_ep = 2.  # set 2m length (long compared to LTX, ensuring we cover all bmvol elements)
				z_endpt = r_ep * np.sin(theta_det_cntr[ith])  # THETA VERTICAL ANGLE
				xy_endpt = r_ep * np.cos(theta_det_cntr[ith])
				# x_endpt, y_endpt = xy_endpt * np.sin(phi_det_cntr[iph]), xy_endpt * np.cos(
				# 	phi_det_cntr[iph])
				x_endpt = x_det + xy_endpt * np.sin(phi_det_cntr[iph])
				y_endpt = y_det + xy_endpt * np.cos(phi_det_cntr[iph])  # PHI MIDPLANE ANGLE
				
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
				
				# det_phi_theta_rmag = [bmvol, phi_tor, [phi,theta,rmag]]
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
				em2d[iph, ith] = em_sum  # todo: divide by solid angle here given dphi&dtheta
		t2 = datetime.datetime.now()
		print(f'... {(t2 - t1).seconds} sec')
		
		print(f'computing flux through detector aperture', end='')
		t1 = datetime.datetime.now()
		if aperture_view is None:
			imax = np.where(em2d == em2d.max())
			aperture_phi0, aperture_theta0 = phi_det_cntr[imax[0][0]], theta_det_cntr[imax[1][0]]
		else:
			aperture_phi0, aperture_theta0 = aperture_view[0], aperture_view[1]
		aperture_radius = 5  # [cm]
		apterture_target_dist = 50  # [cm] (these are guesses for now)
		det_acceptance_angle = np.arctan2(aperture_radius, apterture_target_dist)
		pcirc = np.linspace(0, 2. * np.pi, endpoint=True)
		xcirc, ycirc = aperture_radius * np.cos(pcirc), aperture_radius * np.sin(pcirc)
		phicirc = aperture_phi0 + np.arctan2(xcirc, apterture_target_dist)
		thetacirc = aperture_theta0 + np.arctan2(ycirc, apterture_target_dist)
		
		nphaper, nthaper = 25, 26
		phi_detaper_bins = np.linspace(aperture_phi0 - det_acceptance_angle, aperture_phi0 + det_acceptance_angle,
		                               endpoint=True, num=nphaper + 1)
		phi_detaper_cntr = phi_detaper_bins[0:-1] + (phi_detaper_bins[1] - phi_detaper_bins[0]) / 2.
		theta_detaper_bins = np.linspace(aperture_theta0 - det_acceptance_angle, aperture_theta0 + det_acceptance_angle,
		                                 endpoint=True, num=nthaper + 1)
		theta_detaper_cntr = theta_detaper_bins[0:-1] + (theta_detaper_bins[1] - theta_detaper_bins[0]) / 2.
		em2d_aper = np.zeros((nphaper, nthaper, len(fbm.energy)))  # keep energy distribution here?
		# todo: SHOULD be able to find all elements of det_phi_theta_rmag that fall w/i aperture window in 1 step- can we use rmag to compute the rest? Maybe not if we need to interp psi values and such
		for iph in np.arange(nphaper):  # emissiv_det = [bmvol, iphi, energy]
			for ith in np.arange(nthaper):
				rcheck = np.sqrt(
					(apterture_target_dist * np.tan(phi_detaper_cntr[iph] - aperture_phi0)) ** 2 + (
							apterture_target_dist * np.tan(theta_detaper_cntr[ith] - aperture_theta0)) ** 2)
				if rcheck > aperture_radius:  # point is outside aperture viewing, omit
					em2d_aper[iph, ith, :] = 0.
				else:
					# define path from det outward given angles phi, theta
					# get pts every cm or so along path from det out 2 meters
					r_ep = 2.  # set 2m length (long compared to LTX, ensuring we cover all bmvol elements)
					z_endpt = r_ep * np.sin(theta_detaper_cntr[ith])  # THETA VERTICAL ANGLE
					xy_endpt = r_ep * np.cos(theta_detaper_cntr[ith])
					x_endpt = x_det + xy_endpt * np.sin(phi_detaper_cntr[iph])
					y_endpt = y_det + xy_endpt * np.cos(phi_detaper_cntr[iph])
					
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
					for ipt in np.arange(len(psi_path)):
						# ne_path_m3[ipt] = np.interp(psi_path[ipt], plflx, ne_m3)
						ni_path_m3[ipt] = np.interp(psi_path[ipt], plflx, ni_m3)
					# te_path_ev[ipt] = np.interp(psi_path[ipt], plflx, te_ev)
					# sigv_ei = tabulated_eimpact_ionization(fbm.energy, te_path)
					frac_ioniz = np.outer(ni_path_m3,
					                      (sigv_cx_cm3_s * 1.e-6) / v_ion_m_s * dx)  # vs position and energy
					frac_ioniz = np.nan_to_num(frac_ioniz)  # sets nan to 0.
					frac_trans = 1. - frac_ioniz
					frac_trans = np.cumprod(frac_trans, axis=0)  # cumulative product along position vector v. energy
					
					iwin = np.where(
						(det_phi_theta_rmag[:, :, 0] >= phi_detaper_bins[iph]) & \
						(det_phi_theta_rmag[:, :, 0] < phi_detaper_bins[iph + 1]) & \
						(det_phi_theta_rmag[:, :, 1] >= theta_detaper_bins[ith]) & \
						(det_phi_theta_rmag[:, :, 1] < theta_detaper_bins[ith + 1]))
					
					em_sum = np.zeros_like(fbm.energy)
					for ibmvol in np.arange(len(iwin[0])):
						emiss_v_energy = emissiv_det[iwin[0][ibmvol], iwin[1][ibmvol], :]
						r_to_bmvol = det_phi_theta_rmag[iwin[0][ibmvol], iwin[1][ibmvol], 2]
						iclos = np.where(abs(pts_along - r_to_bmvol) == min(abs(pts_along - r_to_bmvol)))[0][0]
						em_sum += frac_trans[iclos, :] * emiss_v_energy  # keep as array vs energy
					em2d_aper[iph, ith,
					:] = em_sum  # array of #/sec/m^2 vs energy  # todo: need to multiply by solid angle
		em_v_en = np.nansum(em2d_aper, axis=0)
		em_v_en = np.nansum(em_v_en, axis=0)
		t2 = datetime.datetime.now()
		print(f'... {(t2 - t1).seconds} sec')
		
		# create optimal detector alignment array
		ll = np.linspace(0, 1.5)  # pts away from detector
		phi_combined = phi_det2machinecenter - aperture_phi0
		xl, yl = x_det - ll * np.sin(phi_combined), y_det + ll * np.cos(phi_combined)
		phi_detopt, r_detopt = np.arctan2(-xl, yl), np.sqrt(xl ** 2 + yl ** 2)
		
		em2d = np.ma.masked_where(em2d <= 0, em2d)
		
		nr = 50
		rr_bins = np.linspace(0, .7, endpoint=True, num=nr + 1)
		rr = np.linspace(0, .7, endpoint=False, num=nr) + .7 / 50 / 2.
		r2, p2 = np.meshgrid(rr, phi_bin_limits + dphi)
		
		# plot of 3d Halo Neutral Density
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
		
		# top down view of emissivity reaching detector
		em_top = np.zeros((nphi_tor + 1, len(rr)))
		for iph in np.arange(nphi_tor):
			for ir in np.arange(len(rr)):
				em_top[iph, ir] = np.sum(emissiv_det[:, iph, 0][np.where(
					(fbm.r2d * 1.e-2 > rr_bins[ir]) & (fbm.r2d * 1.e-2 < rr_bins[ir + 1])
				)])
				if iph == 0:
					em_top[nphi_tor, ir] = em_top[iph, ir]
		phi_det = np.arctan2(-x_det, y_det)
		em_top = np.ma.masked_where(em_top <= 0, em_top)
		
		savdat = {'fbm': fbm, 'em_v_en': em_v_en, 'p2': p2, 'r2': r2, 'em_top': em_top, 'phi_detopt': phi_detopt,
		          'r_detopt': r_detopt, 'phi_det': phi_det, 'r_det': r_det, 'phi_bb': phi_bb, 'r_bb': r_bb,
		          'n0tot_mp': n0tot_mp, 'phi_det_cntr': phi_det_cntr, 'theta_det_cntr': theta_det_cntr, 'em2d': em2d,
		          'phicirc': phicirc, 'thetacirc': thetacirc}
		pickle.dump(savdat, open(f'{savfn}', 'wb'))
		print(f'saved set data to pkl file {savfn}')
	else:
		savdat = pickle.load(open(f'{savfn}', 'rb'))
		fbm, em_v_en, p2, r2, em_top = savdat['fbm'], savdat['em_v_en'], savdat['p2'], savdat['r2'], savdat['em_top']
		phi_detopt, r_detopt, phi_det = savdat['phi_detopt'], savdat['r_detopt'], savdat['phi_det']
		r_det, phi_bb, r_bb, n0tot_mp = savdat['r_det'], savdat['phi_bb'], savdat['r_bb'], savdat['n0tot_mp']
		phi_det_cntr, theta_det_cntr, em2d = savdat['phi_det_cntr'], savdat['theta_det_cntr'], savdat['em2d']
		phicirc, thetacirc = savdat['phicirc'], savdat['thetacirc']
	
	print('plotting')
	fs = 10
	fig = plt.figure()
	ax1 = fig.add_subplot(221)
	ax2 = fig.add_subplot(222, projection='polar')
	ax3 = fig.add_subplot(224, projection='polar')
	ax4 = fig.add_subplot(223)
	
	# fast-ion distribution plot
	fifig, fiax = plt.subplots(nrows=2)
	fbm.plot(axes=fiax)
	fifig.suptitle('Fast Ion distribution (TRANSP)')
	
	# detector view
	# fig, ax1 = plt.subplots()
	cb1 = ax1.contourf(phi_det_cntr, theta_det_cntr, em2d.transpose(), locator=ticker.LogLocator())
	ax1.plot(phicirc, thetacirc, 'k--')
	cbar1 = fig.colorbar(cb1, ax=ax1)
	cbar1.set_label('num/sec/m^2', fontsize=fs)  # todo: correct units to include per solid angle
	ax1.axvline(0, c='k', ls='--', linewidth=1)
	ax1.axhline(0, c='k', ls='--', linewidth=1)
	# ax.spines['left'].set_position('zero')
	ax1.spines['right'].set_color('none')
	# ax.spines['bottom'].set_position('zero')
	ax1.spines['top'].set_color('none')
	ax1.set_xlabel('$\\phi$ (midplane angle)', fontsize=fs)
	ax1.set_ylabel('$\\theta$ (vertical angle)', fontsize=fs)
	ax1.axis('equal')
	
	# top down view of emissivity reaching detector
	# fig2, ax2 = plt.subplots(subplot_kw=dict(projection='polar'))
	cb2 = ax2.contourf(p2, r2, em_top, locator=ticker.LogLocator())
	if tangential:
		ax2.plot(phi_detopt, r_detopt, 'r--')
	cbar2 = fig.colorbar(cb2, ax=ax2)
	cbar2.set_label('num/sec/m^2',
	                fontsize=fs)  # units correct- no solid angle in this plot, just plotting rate/m^2 onto detector
	ax2.plot(phi_det, r_det, 'ro')
	ax2.annotate('detector', (phi_det, r_det), fontsize=fs)
	ax2.set_theta_zero_location('N')
	
	# fig3, ax3 = plt.subplots(subplot_kw=dict(projection='polar'))
	cb3 = ax3.contourf(phi_bb, r_bb, n0tot_mp, levels=25)
	if tangential:
		ax3.plot(phi_detopt, r_detopt, 'r--')
	ax3.plot(phi_det, r_det, 'ro')
	cbar3 = fig.colorbar(cb3, ax=ax3)
	cbar3.set_label('#/cm^3\n(summed over vertical coord)', fontsize=fs)
	rlim = ax3.get_ylim()
	ax3.set_rlim((0, rlim[1]))
	ax3.set_theta_zero_location('N')
	ax3.set_title('beam neutral halo', fontsize=fs)
	
	# fig4, ax4 = plt.subplots()
	ax4.bar(fbm.energy / 1000., em_v_en, width=.9)
	ax4.set_xlim(right=20)
	ax4.set_xlabel('incident ion energy (keV)', fontsize=fs)
	ax4.set_ylabel('flux onto detector (#/sec)', fontsize=fs)
	
	for ax in [ax1, ax2, ax3, ax4]:
		ax.tick_params(labelsize=fs)
	
	plt.tight_layout()


def set_analysis(shot):
	'''
	designed specifically for transp runs 106536R02,3,4 and 105795R02,3,4
	2,3,4 vary by tangency radius 2=24cm, 3=21.3, 4=33
	'''
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	fig1, ((ax11, ax12), (ax13, ax14)) = plt.subplots(nrows=2, ncols=2, sharex='col', figsize=[10, 5])
	rtans = [21.3, 24, 33]
	npa = NPA()  # attempt to pull in NPA instrument function- not yet relevant data
	for iset, rtan in enumerate([3, 2, 4]):
		ax11.plot([np.nan], [np.nan], label=f'{rtans[iset]}', c=clrs[iset])
		pkl_fn = f'{direc}t{shot}/{shot}R{rtan:02}.pkl'
		dat = pickle.load(open(pkl_fn, 'rb'))
		fbm = dat['fbm']
		fbeam = np.zeros((len(fbm.energy), len(fbm.pitch)))
		for ie in np.arange(len(fbm.energy)):
			for ip in np.arange(len(fbm.pitch)):
				fbeam[ie, ip] = np.mean(fbm.fbm[:, ip, ie])
		fbeam[fbeam <= 0] = 1  # set positive minimum to avoid issues while plotting on log scale
		# levels = np.logspace(0., np.log10(np.max(fbeam)), num=25)
		fbeam /= np.nanmax(fbeam)
		ch_energies, ch_curr = npa.compute_channel_current(fbm.energy / 1000., dat['em_v_en'])  # [A], [kV]
		
		levels = [.2, .4, .6, .8, .9]
		cont = ax11.contour(fbm.energy / 1000., fbm.pitch, fbeam.transpose(), locator=ticker.LogLocator(),
		                    levels=levels, colors=clrs[iset])
		ax11.clabel(cont, fontsize=9, inline=True)
		ax13.plot(fbm.energy / 1000., fbm.n_tot_arr, c=clrs[iset])
		ax12.bar(fbm.energy / 1000., dat['em_v_en'], width=.9 * (fbm.energy[1] - fbm.energy[0]) * 1.e-3,
		         edgecolor=clrs[iset], fill=False)
		ax14.plot(ch_energies, ch_curr, 'o-', c=clrs[iset])
		if iset == 2:
			chtext = ['ch# 0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
			for ich in range(len(chtext)):
				ax14.annotate(chtext[ich], (ch_energies[ich], ch_curr[ich]))
	
	ax11.set_ylim((0, 1))
	ax11.set_ylabel('Pitch')
	for ax in [ax13, ax14]:
		ax.set_xlabel('energy (keV)')
		ax.set_xlim((0, 20))
	ax12.set_ylabel('NPA flux (#/sec)')
	ax13.set_ylabel('# per energy bin')
	ax14.set_ylabel('channel current (A)')
	ax11.legend(title='rtan (cm)')
	fig1.suptitle(f'{shot}')
	plt.tight_layout()


if __name__ == '__main__':
	direc = 'Z:/transp/'
	shot = 105795
	# shot = 106536
	
	create_npa_data = 1
	# for shot in [105795, 106536]:
	# 	for rtan in [2, 3, 4]:
	if create_npa_data:
		for shot in [106536]:
			# for (rtanfn, rtan) in zip([2, 3, 4], [24., 21.3, 33]):
			for (rtanfn, rtan) in zip([2], [24.]):
				fbm_fn = f'{direc}t{shot}/{shot}R{rtanfn:02}_fi_5.cdf'
				n0_fn = f'{direc}t{shot}/{shot}R{rtanfn:02}_boxn0_5.cdf'
				eq_fn = f'{direc}t{shot}/{shot}R{rtanfn:02}_05.eqdsk'
				nenite_fn = f'{direc}t{shot}/{shot}R{rtanfn:02}_05.nenite'
				synthetic_npa(fbm_fn, n0_fn, eq_fn, nenite_fn, savfn=f'{direc}t{shot}/{shot}R{rtanfn:02}radial',
				              beamtan=rtan)
	
	do_set_analysis = 0
	if do_set_analysis:
		set_analysis(shot)
	
	plt.show()
