"""
author: wcapecch
These classes read in data from .CDF files
The CDF files are created by the user after a TRANSP run completes
Details on generating the CDF files can be found on the LTXb Wiki
"""
import datetime

from matplotlib import ticker
import numpy as np
import matplotlib.pyplot as plt
import MDSplus
from scipy.io import netcdf
from helpful_stuff import read_transp_cdf, closest


class FBM:
	def __init__(self, cdf_file, label=None):
		tstrt = datetime.datetime.now()
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
		tstp = datetime.datetime.now()
		print(f'FBM cdf load time: {(tstp-tstrt).seconds} sec')
		
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
		print(f'Reading in {cdf_file}... be patient')
		tstrt = datetime.datetime.now()
		cdf = netcdf.netcdf_file(cdf_file, 'r', mmap=False)
		vars = cdf.variables
		self.vars = vars
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
		self.boxn0 = vars['BOXN0'].data  # [NLBOX, NYBOX, NXBOX, NBBOX] (reverse of in desc)
		self.boxn0_desc = vars['BOXN0'].BOXN0
		self.boxn0h0 = vars['BOXN0H0'].data
		self.boxn0h0_desc = vars['BOXN0H0'].BOXN0H0
		self.boxn0hh = vars['BOXN0HH'].data
		self.boxn0hh_desc = vars['BOXN0HH'].BOXN0HH
		self.dx = self.xbox[1] - self.xbox[0]  # cm
		self.dy = self.ybox[1] - self.ybox[0]  # cm
		self.dl = self.lbox[1] - self.lbox[0]  # cm
		self.dv = self.dx * self.dy * self.dl  # cm^3
		self.n0_vs_lbox = np.sum(self.boxn0[:, :, :, 0], axis=(1, 2))  # sum over x, y  #/cm^3 vs lbox
		self.num_vs_lbox = self.n0_vs_lbox * self.dv  # num vs lbox
		cdf.close()
		print('N0_tot = {:.4e}'.format(np.sum(self.num_vs_lbox)))
		tstp = datetime.datetime.now()
		print(f'Halo neutral CDF load time {(tstp - tstrt).seconds} sec')
	
	def plot(self, xslice=None, yslice=None, lslice=None, label=None, nbeam=1):
		nplots = 1
		if xslice is not None:
			nplots += 1
		if yslice is not None:
			nplots += 1
		if lslice is not None:
			nplots += 1
		fig = plt.figure(self.cdf_file)
		ax = fig.add_subplot(1, nplots, 1)
		ax.plot(self.lbox, self.n0_vs_lbox)
		ax.set_xlabel('LBOX (cm)')
		ax.set_ylabel('#/cm^3')
		iplot = 2
		if xslice is not None:
			val, ix = closest(self.xbox, xslice)
			xax = fig.add_subplot(1, nplots, iplot, projection='3d')
			iplot += 1
			yy, ll = np.meshgrid(self.ybox, self.lbox)
			xax.plot_surface(yy, ll, self.boxn0[:, :, ix, nbeam - 1], cmap='viridis')
			xax.set_zlabel('n0 (#/cm^3)')
			xax.set_xlabel('YBOX (cm)')
			xax.set_ylabel('LBOX (cm)')
			xax.set_title(f'Slice at XBOX={val}cm')
		if yslice is not None:
			val, iy = closest(self.ybox, yslice)
			yax = fig.add_subplot(1, nplots, iplot, projection='3d')
			iplot += 1
			xx, ll = np.meshgrid(self.xbox, self.lbox)
			yax.plot_surface(xx, ll, self.boxn0[:, iy, :, nbeam - 1], cmap='viridis')
			yax.set_zlabel('n0 (#/cm^3)')
			yax.set_xlabel('XBOX (cm)')
			yax.set_ylabel('LBOX (cm)')
			yax.set_title(f'Slice at YBOX={val}cm')
		if lslice is not None:
			val, il = closest(self.lbox, lslice)
			lax = fig.add_subplot(1, nplots, iplot, projection='3d')
			xx, yy = np.meshgrid(self.xbox, self.ybox)
			lax.plot_surface(xx, yy, self.boxn0[il, :, :, nbeam - 1], cmap='viridis')
			lax.set_zlabel('n0 (#/cm^3)')
			lax.set_xlabel('XBOX (cm)')
			lax.set_ylabel('YBOX (cm)')
			lax.set_title(f'Slice at LBOX={val}cm')


class SimpleCDF:
	"""Structured in a compatible way to SimpleSignal, we look for TRANSP data in a CDF file
	"""
	
	def __init__(self, shot, cdf, name):
		local_dir = '//samba/wcapecch/transp_rawdata/'
		
		# Initialize the fields of the simpleSignal class to be the input of
		# the constructor.
		self.shot = shot
		self.nodepath = name
		self.tree = 'local_CDF_file'
		self.data = None
		self.dataunits = None
		self.ndim = None
		self.dim1 = None
		self.dim1units = None
		self.dim2 = None
		self.dim2units = None
		self.dim3 = None
		self.dim3units = None
		self.error = None
		self.name = None
		self.label = None
		self.goodshot = None  # indicates whether or not a shot exists in tree with any data
		
		if name.upper() in cdf.keys():
			dat = cdf[name.upper()]
			self.data = dat.data  # .transpose()  # need transpose here for some reason for good plotting
			self.dataunits = dat.units
			self.ndim = dat.data.ndim
			if self.ndim >= 1:
				dim1 = cdf[dat.dimensions[0]]
				if dim1.data.ndim == 1:
					self.dim1 = dim1.data
					self.dim1units = dim1.units
				else:  # dim1 has multiple dimensions
					# ASSUME FOR NOW IT ALWAYS GOES:: time, data
					islice = np.where(np.array(dim1.dimensions) == dat.dimensions[0])[0][0]
					if islice != 1:
						a = 1
					self.dim1 = dim1.data[0, :]
					self.dim1units = dim1.units
			if self.ndim >= 2:
				# NOTE for some reason coming from MDSPlus we get time on 2nd dimension
				# So here we're flipping to correct for difference between MDSPlus and CDF data pull
				self.dim2 = self.dim1
				self.dim2units = self.dim1units
				dim1 = cdf[dat.dimensions[1]]
				if dim1.data.ndim == 1:
					self.dim1 = dim1.data
					self.dim1units = dim1.units
				else:  # dim1 has multiple dimensions
					# ASSUME FOR NOW IT ALWAYS GOES:: time, data
					islice = np.where(np.array(dim1.dimensions) == dat.dimensions[1])[0][0]
					if islice != 1:
						a = 1
					self.dim1 = dim1.data[0, :]
					self.dim1units = dim1.units
			if self.ndim >= 3:
				a = 1  # what signal is this??
			self.error = None
			self.name = name
			self.label = name
			self.goodshot = True
		else:
			self.goodshot = False
	
	def plot(self, title=None, contour=False, label=None, ax=None, fontsize=None):
		"""Plot the signal vs. time using as much of the stored
		 information as is available to annotate the figure. Set the
		 color of the plot line using, e.g. color='b' for blue."""
		
		if ax is None:
			fig = plt.figure()
		if fontsize is not None:
			plt.rcParams.update({'font.size': fontsize})
		if self.ndim == 1:
			plt.plot(self.dim1, self.data, label=label)
			plt.xlabel('x ({})'.format(self.dim1units))
			plt.ylabel('{} ({})'.format(self.label, self.dataunits))
		elif self.ndim == 2:
			if contour:
				if ax is None:
					ax = plt.axes()
				ax.contourf(self.dim1, self.dim2, self.data)
				plt.title('{} ({})'.format(self.label, self.dataunits))
			else:
				if ax is None:
					ax = plt.axes(projection='3d')
				xx, yy = np.meshgrid(self.dim1, self.dim2)
				ax.plot_surface(xx, yy, self.data, cmap='viridis')
				ax.set_zlabel('{} ({})'.format(self.label, self.dataunits))
			ax.set_xlabel('x ({})'.format(self.dim1units))
			ax.set_ylabel('y ({})'.format(self.dim2units))
			ax.set_title(title)
		elif self.ndim == 3:
			num = len(self.data)  # len here gives length of 1st dimension
			xx, yy = np.meshgrid(self.dim2, self.dim3)  # assume dim1 is time axis
			iplot = [int(ii) for ii in np.linspace(0, len(self.dim1) - 1, num=3)]
			for n, i in enumerate(iplot):
				if contour:
					ax = fig.add_subplot(3, 1, n + 1)
					ax.contourf(self.dim2, self.dim3, self.data[:, :, i],
					            label='{} ({})'.format(self.label, self.dataunits))
					ax.legend()
				else:
					ax = fig.add_subplot(3, 1, n + 1, projection='3d')
					ax.plot_surface(xx, yy, self.data[:, :, i], cmap='viridis')
					ax.set_zlabel('{} ({})'.format(self.label, self.dataunits))
					ax.set_title('t: {}'.format(self.dim1[i]))
				ax.set_xlabel('x ({})'.format(self.dim1units))
				ax.set_ylabel('y ({})'.format(self.dim2units))


if __name__ == '__main__':
	do_fbm_example = 0
	if do_fbm_example:
		# EXAMPLE usage of FBM class
		fn = 'Z:/transp/t111539/111539F04_fi_3.cdf'
		fbm = FBM(fn)
		fbm.plot(label='')
	
	do_halo3d_example = 1
	if do_halo3d_example:
		# EXAMPLE usage of Halo3D class
		fn = 'Z:/transp_rawdata/103617C03_boxn0_2.cdf'
		halo3d = Halo3D(fn)
		halo3d.plot(label='', xslice=0.1, yslice=-5., lslice=200.)
		plt.show()
	
	do_simplecdf_example = 0
	if do_simplecdf_example:
		# EXAMPLE usage of SimpleCDF
		run = 1000022029
		cdf = read_transp_cdf(run)
		if cdf is not None:
			bsorb_h = SimpleCDF(run, cdf, 'BSORB_H')
			sbdepmc_h = SimpleCDF(run, cdf, 'SBDEPMC_H')
			shinethru = SimpleCDF(run, cdf, 'SBSHINE_H')
			sinj = SimpleCDF(run, cdf, 'SINJ')
			bdens = SimpleCDF(run, cdf, 'BDENS2_H')  # used in Bill's method of coupled fraction
