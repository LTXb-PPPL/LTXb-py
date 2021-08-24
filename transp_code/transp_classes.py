"""
author: wcapecch
These classes read in data from .CDF files
The CDF files are created by the user after a TRANSP run completes
Details on generating the CDF files can be found on the LTXb Wiki
"""

from matplotlib import ticker
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt


class FBM:
	def __init__(self, cdf_file, label=None):
		cdf = netcdf.netcdf_file(cdf_file, 'r', mmap=False)
		vars = cdf.variables
		self.cdf_file = cdf_file
		self.label = label
		self.energy = vars['E_H_NBI'].data
		self.energyunits = vars['E_H_NBI'].units
		self.pitch = vars['A_H_NBI'].data
		self.fbm = vars['F_H_NBI'].data
		self.fbmunits = vars['F_H_NBI'].units
		self.r2d = vars['R2D'].data
		self.r2dunits = vars['R2D'].units
		self.z2d = vars['Z2D'].data
		self.z2dunits = vars['Z2D'].units
		self.bmvol = vars['BMVOL'].data
		self.bmvolunits = vars['BMVOL'].units
		self.n_tot, self.n_tot_arr = self.compute_n_tot()  # vs energy
		cdf.close()
	
	def compute_n_tot(self):
		# fbm has indices of [bmvol, pitch, energy]
		n_tot_arr = np.zeros_like(self.energy)
		e_av = 0.
		for i in np.arange(len(self.energy) - 1):
			n = 0.
			de = self.energy[i + 1] - self.energy[i]
			for k in np.arange(len(self.bmvol)):
				n += np.sum(self.fbm[k, :, i] * self.bmvol[k] * de / len(self.pitch))
				e_av += np.sum(self.fbm[k, :, i] * self.bmvol[k] * self.energy[i] * de / len(self.pitch))
			n_tot_arr[i] = n
		e_av /= np.sum(n_tot_arr)
		print('N_tot = {:.4e}'.format(np.sum(n_tot_arr)))
		print('E_av (eV) = {:.4e}'.format(e_av))
		return np.sum(n_tot_arr), n_tot_arr
	
	def compute_num(self, r_cm=None, z_cm=None, vs_pitch=True):
		# fbm has indices of [bmvol, pitch, energy]
		if r_cm is not None and z_cm is not None:
			wzones = np.where((self.r2d - r_cm) ** 2 + (
					self.z2d - z_cm) ** 2 <= 2. ** 2.)  # we're copying deyong here who used r=10.
		else:
			wzones = np.arange(len(self.r2d))  # all the indices
		if vs_pitch:
			xx = self.pitch
			n_tot_arr = np.zeros_like(self.pitch)
			de = np.mean([self.energy[e + 1] - self.energy[e] for e in np.arange(len(self.energy) - 1)])
			e_av = 0.
			for p in np.arange(len(self.pitch) - 1):
				n = 0.
				for k in wzones[0]:
					n += np.sum(self.fbm[k, p, :] * self.bmvol[k] * de / len(self.energy))
					e_av += np.sum(self.fbm[k, p, :] * self.bmvol[k] * self.energy * de / len(self.energy))
				n_tot_arr[p] = n
			e_av /= np.sum(n_tot_arr)
			print('N_tot = {:.4e}'.format(np.sum(n_tot_arr)))
			print('E_av (eV) = {:.4e}'.format(e_av))
		else:
			xx = self.energy
			n_tot_arr = np.zeros_like(self.energy)
			e_av = 0.
			for i in np.arange(len(self.energy) - 1):
				n = 0.
				de = self.energy[i + 1] - self.energy[i]
				for k in np.arange(len(self.bmvol)):
					n += np.sum(self.fbm[k, :, i] * self.bmvol[k] * de / len(self.pitch))
					e_av += np.sum(self.fbm[k, :, i] * self.bmvol[k] * self.energy[i] * de / len(self.pitch))
				n_tot_arr[i] = n
			e_av /= np.sum(n_tot_arr)
			print('N_tot = {:.4e}'.format(np.sum(n_tot_arr)))
			print('E_av (eV) = {:.4e}'.format(e_av))
		return xx, n_tot_arr
	
	def plot(self, r_cm=None, z_cm=None, label=None, axes=None):
		if axes is None:
			fig = plt.figure(self.cdf_file)
			# plt.suptitle(label)
			ax1, ax2 = fig.add_subplot(211), fig.add_subplot(212)
			subplotting = False
		else:
			ax1, ax2 = axes[0], axes[1]
			subplotting = True
		ax1.set_title(label)
		if r_cm is not None and z_cm is not None:
			wzones = np.where(
				(self.r2d - r_cm) ** 2 + (self.z2d - z_cm) ** 2 <= 2. ** 2.)  # we're copying deyong here who used r=10.
		else:
			wzones = np.arange(len(self.r2d))  # all the indices
		rplot = np.mean(self.r2d[wzones])
		zplot = np.mean(self.z2d[wzones])
		print('Plot the fast-ion density distribution at (r, z) = ({}, {})'.format(rplot, zplot))
		fbeam = np.zeros((len(self.energy), len(self.pitch)))
		# rel_num = np.zeros(len(self.energy))
		for ie in np.arange(len(self.energy)):
			for ip in np.arange(len(self.pitch)):
				fbeam[ie, ip] = np.mean(self.fbm[wzones, ip, ie])
		# rel_num[ie] = np.mean(self.fbm[wzones, :, ie])
		fbeam[fbeam <= 0] = 1  # set positive minimum to avoid issues while plotting on log scale
		levels = np.logspace(0., np.log10(np.max(fbeam)), num=25)
		cont = ax1.contourf(self.energy, self.pitch, fbeam.transpose(), locator=ticker.LogLocator(), levels=levels)
		# cbar = fig.colorbar(cont, ax=ax1)
		ax2.plot(self.energy, self.n_tot_arr, 'o-')
		if not subplotting:
			ax1.set_ylabel('Pitch')
			ax2.set_xlabel('energy (eV)')
			ax2.set_ylabel('relative number')
		ax1.set_xlim((0, 30.e3))
		ax2.set_xlim((0, 30.e3))


class Halo3D:
	def __init__(self, cdf_file, label=None):
		cdf = netcdf.netcdf_file(cdf_file, 'r', mmap=False)
		vars = cdf.variables
		self.cdf_file = cdf_file
		self.label = label
		self.numavg = vars['NUMAVG'].data
		self.numavg_desc = vars['NUMAVG'].NUMAVG
		self.avgtim = vars['AVGTIM'].data
		self.avgtim_desc = vars['AVGTIM'].AVGTIM
		self.mthdavg = vars['MTHDAVG'].data
		self.mthdavg_desc = vars['MTHDAVG'].MTHDAVG
		self.avgsamp = vars['AVGSAMP'].data
		self.avgsamp_desc = vars['AVGSAMP'].AVGSAMP
		self.abeama = vars['ABEAMA'].data
		self.abeama_desc = vars['ABEAMA'].ABEAMA
		self.zbeama = vars['ZBEAMA'].data
		self.zbeama_desc = vars['ZBEAMA'].ZBEAMA
		self.einja = vars['EINJA'].data
		self.einja_desc = vars['EINJA'].EINJA
		self.xbox = vars['XBOX'].data
		self.xbox_desc = vars['XBOX'].XBOX
		self.ybox = vars['YBOX'].data
		self.ybox_desc = vars['YBOX'].YBOX
		self.lbox = vars['LBOX'].data
		self.lbox_desc = vars['LBOX'].LBOX
		self.boxn0 = vars['BOXN0'].data
		self.boxn0_desc = vars['BOXN0'].BOXN0
		self.boxn0h0 = vars['BOXN0H0'].data
		self.boxn0h0_desc = vars['BOXN0H0'].BOXN0H0
		self.boxn0hh = vars['BOXN0HH'].data
		self.boxn0hh_desc = vars['BOXN0HH'].BOXN0HH
		self.dx = self.xbox[1] - self.xbox[0]  # cm
		self.dy = self.ybox[1] - self.ybox[1]  # cm
		self.dl = self.lbox[1] - self.lbox[0]  # cm
		self.dv = self.dx * self.dy * self.dl  # cm^3
		self.n0_vs_lbox = np.sum(self.boxn0[:, :, :, 0], axis=(1, 2))  # sum over x, y  #/cm^3 vs lbox
		self.num_vs_lbox = self.n0_vs_lbox * self.dv  # num vs lbox
		cdf.close()
		print('N0_tot = {:.4e}'.format(np.sum(self.num_vs_lbox)))
	
	def plot(self, label=None):
		fig = plt.figure(self.cdf_file)
		plt.plot(self.lbox, self.n0_vs_lbox)
		plt.xlabel('LBOX (cm)')
		plt.ylabel('#/cm^3')


if __name__ == '__main__':
	# EXAMPLE usage of FBM class
	fn = 'Z:/wcapecch/transp/t111539/111539F04_fi_3.cdf'
	fbm = FBM(fn)
	fbm.plot(label='')
	
	# EXAMPLE usage of Halo3D class
	fn = 'Z:/wcapecch/transp_rawdata/103617C02_boxn0_5.cdf'
	halo3d = Halo3D(fn)
	halo3d.plot(label='')
