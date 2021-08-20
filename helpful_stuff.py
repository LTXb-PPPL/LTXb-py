try:
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt
except:
	pass
import numpy as np
import MDSplus
import re
import datetime
from bills_LTX_MDSplus_toolbox import get_tree_conn
import glob
import os


def get_current_ltx_shot():
	# latest shot in 'ltx_b' tree
	shotdone = open('Y:\\maindata\\shotdone', 'r')
	shot = int(shotdone.read())
	
	return shot


def get_current_nbi_shot(ltx_shotnums=False):
	# this is used when whether or not last shot was paired with LTX discharge is unknown
	# so this will give most recent nbi discharge whether it was NBI-only or NBI+LTX
	
	# old method: this gets largest shot number from tree, which causes issues since LTX shotnums are smaller than NBI ones
	# latest shot in 'ltx_nbi' tree
	# tc = get_tree_conn(0, 'ltx_nbi')  # open tree connection
	# time = tc.get('.metadata.timestamp')  # grab timestamp of last shot
	# dt = datetime.datetime.strptime(time.split(' ')[0], '%m/%d/%Y')
	# fldr = dt.strftime('%y%m%d')
	
	date_dirs = [d for d in os.listdir('Y:\\NBI\\data\\') if d.isnumeric()]
	date_dirs.sort()
	for fldr in date_dirs[::-1]:  # go from most recent (format YYMMDD)
		path = f'Y:\\NBI\\data\\{fldr}\\'
		files = glob.glob(path + 'nbi_adc*.tdms')
		shotnums = [f.split('.tdms')[0] for f in files]
		shotnums = [int(sh.split('nbi_adc')[1]) for sh in shotnums]
		if len(files) == 0:
			print('\nNO FILES FOUND: open the Y: drive to connect for some reason')
			return np.nan
		else:
			if ltx_shotnums:
				shots = [sh for sh in shotnums if sh < 200000]
				if len(shots) > 0:
					return max(shots)
			else:
				shots = [sh for sh in shotnums if sh > 200000]
				if len(shots) > 0:
					return max(shots)


def read_nenite(nenite_fn):
	# These files are created in dump_denstemp_from_cdf.pro on NoMachine
	# :param nenite_fn:
	# :return:
	# 	xb : sqrt(phi/philim) Phi is the toroidal flux enclosed within a surface, and Philim the toroidal flux enclosed within the boundary
	# 	plflx : Poloidal Flux [wb/rad]
	# 	ne  : electron density [n/cm^3] interpolated from \x onto \xb
	# 	ni  : ion density [n/cm^3] interpolated from \x onto \xb
	# 	te  : electron temp [eV] interpolated from \x onto \xb
	f = open(nenite_fn, 'r')
	form = r'-?\d+.\d+E[+-]\d+'
	for i in np.arange(6):  # go through header
		f.readline()
	nn = int(re.findall(r'\d+', f.readline())[0])
	nrows = np.ceil(nn / 6.)  # num rows to scan fpol, etc
	xb, plflx, ne, ni, te = [], [], [], [], []
	for i in np.arange(nrows) + 1:
		xb.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		plflx.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		ne.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		ni.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		te.extend([float(n) for n in re.findall(form, f.readline())])
	return np.array(xb), np.array(plflx), np.array(ne), np.array(ni), np.array(te)


def read_eqdsk(eqdsk, plot=False):
	f = open(eqdsk, 'r')
	form = r'\S\d+.\d+E[+-]\d+'
	_, nxeqd, nyeqd = [int(dig) for dig in re.findall(r'\d+', f.readline())]
	xdimeqd, ydimeqd, r0, redeqd, ymideqd = [float(n) for n in re.findall(form, f.readline())]
	xma, yma, psimag, psilim, beqd = [float(n) for n in re.findall(form, f.readline())]
	toteqd, psimx1, psimx2, xax1, xax2 = [float(n) for n in re.findall(form, f.readline())]
	zax1, zax2, psisep, xsep, ysep = [float(n) for n in re.findall(form, f.readline())]
	nrows = np.ceil(nxeqd / 5.)  # num rows to scan fpol, etc
	fpol, pres, ffpeqd, ppeqd, psixy, q = [], [], [], [], [], []
	for i in np.arange(nrows) + 1:
		fpol.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		pres.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		ffpeqd.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		ppeqd.extend([float(n) for n in re.findall(form, f.readline())])
	for j in np.arange(nyeqd) + 1:
		psirow = []
		for i in np.arange(nrows) + 1:
			psirow.extend([float(n) for n in re.findall(form, f.readline())])
		psixy.append(psirow)
	for i in np.arange(nrows) + 1:
		q.extend([float(n) for n in re.findall(form, f.readline())])
	nrlimit, nrves = [int(dig) for dig in re.findall(r'\d+', f.readline())]
	nrows_lim = np.ceil(2 * nrlimit / 5.)
	nrows_ves = np.ceil(2 * nrves / 5)
	rzlim, rzves = [], []
	for i in np.arange(nrows_lim) + 1:
		rzlim.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows_ves) + 1:
		rzves.extend([float(n) for n in re.findall(form, f.readline())])
	f.close()
	rlimit, zlimit = [rzlim[2 * i] for i in np.arange(nrlimit)], [rzlim[2 * i + 1] for i in np.arange(nrlimit)]
	rves, zves = [rzves[2 * i] for i in np.arange(nrves)], [rzves[2 * i + 1] for i in np.arange(nrves)]
	
	fpol, pres, ffpeqd, ppeqd, psixy, q = np.array(fpol), np.array(pres), np.array(ffpeqd), np.array(
		ppeqd), np.array(psixy), np.array(q)
	rlimit, zlimit, rves, zves = np.array(rlimit), np.array(zlimit), np.array(rves), np.array(zves)
	
	dR = xdimeqd / (nxeqd - 1.)
	dZ = ydimeqd / (nyeqd - 1.)
	x = np.arange(redeqd, redeqd + xdimeqd + dR, dR)
	y = np.arange(ymideqd - ydimeqd / 2, ymideqd + ydimeqd / 2 + dZ, dZ)
	
	[x_xy, y_xy] = np.meshgrid(x, y)
	
	psiout = np.linspace(psimag, psilim, nxeqd)
	real_psilim = (psilim - .01 * psimag) / .99
	psirzvec = psixy.reshape(nxeqd * nyeqd)
	q_xy = np.interp(psirzvec, psiout, q)
	BtR = np.interp(psirzvec, psiout, fpol)
	BtR = BtR.reshape([nxeqd, nyeqd])
	# bphi_xy = BtR / np.transpose(x_xy)  # WJC commented out 6/23/21- why the transpose?? looks wrong
	bphi_xy = BtR / x_xy
	bphi_mp = fpol / x
	
	rvec_hd, zvec_hd = np.linspace(x[0], x[-1], num=1000.), np.linspace(y[0], y[-1], num=1000)
	r2dhd, z2dhd = np.meshgrid(rvec_hd, zvec_hd)
	nz0 = np.where(abs(y) == min(abs(y)))[0][0]
	br_xy, bz_xy = np.zeros_like(psixy), np.zeros_like(psixy)
	for i in np.arange(len(x)):
		bz_xy[i, :] = np.gradient(psixy[i, :], x) / x
		br_xy[:, i] = -np.gradient(psixy[:, i], y - y[0]) / x[i]
	bp_xy = np.sqrt(br_xy ** 2 + bz_xy ** 2)
	bp_mp, br_mp, bz_mp = bp_xy[nz0, :], br_xy[nz0, :], bz_xy[nz0, :]
	bave_xy = np.sqrt(bp_xy ** 2 + bphi_xy ** 2)
	bave_mp = bave_xy[nz0, :]
	
	if plot:
		plt.plot(x, br_xy[nz0, :], label='br_xy')
		plt.plot(x, bz_xy[nz0, :], label='bz_xy')
		plt.plot(x, bphi_mp, label='bt')
		plt.legend()
		plt.show()
	# BE MORE CAREFUL ABOUT LIMITER HERE
	# for i in np.arange(len(x)):
	# 	for j in np.arange(len(y)):
	# 		if np.sqrt((x[i] - .4) ** 2 + y[j] ** 2) > .3:
	# 			br_xy[j, i], bz_xy[j, i] = np.nan, np.nan
	if x[np.where(bz_mp == min(bz_mp))][0] < x[np.where(bz_mp == max(bz_mp))][0]:
		ip_is_cw = True
	else:
		ip_is_cw = False
	if bphi_mp[0] < 0:
		bt_is_cw = True
	else:
		bt_is_cw = False
	
	output = {'nxeqd': nxeqd, 'nyeqd': nyeqd, 'xdimeqd': xdimeqd, 'ydimeqd': ydimeqd, 'r0': r0, 'redeqd': redeqd,
	          'ymideqd': ymideqd, 'xma': xma, 'yma': yma, 'psimag': psimag, 'psilim': psilim, 'beqd': beqd,
	          'toteqd': toteqd, 'psimx1': psimx1, 'psimx2': psimx2, 'xax1': xax1, 'xax2': xax2, 'zax1': zax1,
	          'zax2': zax2, 'psisep': psisep, 'xsep': xsep, 'ysep': ysep, 'fpol': fpol, 'pres': pres, 'ffpeqd': ffpeqd,
	          'ppeqd': ppeqd, 'q': q, 'nrlimit': nrlimit, 'nrves': nrves, 'rlimit': rlimit, 'zlimit': zlimit,
	          'rves': rves, 'zves': zves, 'x': x, 'y': y, 'psixy': psixy, 'real_psilim': real_psilim, 'q_xy': q_xy,
	          'bphi_xy': bphi_xy, 'bp_xy': bp_xy, 'x_xy': x_xy, 'y_xy': y_xy, 'br_xy': br_xy, 'bz_xy': bz_xy,
	          'bp_mp': bp_mp, 'bphi_mp': bphi_mp, 'bave_mp': bave_mp, 'bave_xy': bave_xy, 'br_mp': br_mp,
	          'bz_mp': bz_mp, 'ip_is_cw': ip_is_cw, 'bt_is_cw': bt_is_cw}
	
	# fpol_xy: fpol_xy, rminor_xy: rminor_xy, norm_psixy: norm_psixy, plasma: new_plasma, dpsi: dpsi,
	# psiaxis: psi_axis}
	return output


def ltx_limiter():
	'''
	myFile = fopen('ltx_limiter.txt','r');
	myNum = fscanf(myFile,'%d',1);
	rlimiter = fscanf(myFile,'%f',myNum);
	zlimiter = fscanf(myFile,'%f',myNum);
	nlim = numel(rlimiter);
	fclose(myFile);
	'''
	lim_fn = '//samba/wcapecch/python/ltx_limiter.txt'
	f = open(lim_fn, 'r')
	form = r'-?\d+.\d+E[+-]\d+'
	nn = int(re.findall(r'\d+', f.readline())[0])
	nrows = np.ceil(nn / 6.)  # num rows to scan fpol, etc
	rlimiter, zlimiter = [], []
	for i in np.arange(nrows) + 1:
		rlimiter.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		zlimiter.extend([float(n) for n in re.findall(form, f.readline())])
	rlimiter, zlimiter = rlimiter[0:-1], zlimiter[0:-1]  # get rid of double pt (doing it better below)
	
	rlimiter = np.roll(rlimiter, -2)
	zlimiter = np.roll(zlimiter, -2)
	r0, z0, r1, z1 = rlimiter[0], zlimiter[0], rlimiter[-1], zlimiter[-1]
	# create double overlap (useful for interpolation arrays generated below)
	rlimiter, zlimiter = np.insert(rlimiter, 0, r1), np.insert(zlimiter, 0, z1)
	rlimiter, zlimiter = np.append(rlimiter, r0), np.append(zlimiter, z0)
	
	rmag = 0.4  # doesn't need to be exact, just using interior pt. to find pts outside limiter
	rminor_lim = np.sqrt((rlimiter - rmag) ** 2 + zlimiter ** 2)
	theta_lim = np.arctan2(zlimiter, (rlimiter - rmag))
	theta_lim[0] -= 2. * np.pi
	theta_lim[-1] += 2. * np.pi
	
	return rlimiter, zlimiter, rminor_lim, theta_lim


def color_cycle_demo():
	# i always forget how to use the color cycle
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	for i, col in enumerate(clrs):
		plt.plot(np.arange(10), np.arange(10) ** i / 9. ** i, col)
	cmap = plt.get_cmap('viridis')
	for i, col in enumerate(cmap.colors):
		plt.plot(-np.arange(10), np.ones(10) * i / len(cmap.colors), col)
	# plt.plot(np.arange(10), clrs[0])
	# plt.plot(np.arange(10) * 1.5, clrs[1])
	a = 1


def plot_deposition_fractions(x, runs, ax, beamon=None, beamduration=None, overplot=0, label='', plot=True):
	if beamon is None:
		beamon = np.zeros_like(runs)
	if beamduration is None:
		beamduration = np.ones_like(runs) * 30.e-6
	coupl_arr = []  # num coupled @ beamoff
	pr_loss_arr = []  # num prompt loss @ beamoff
	num_inj_arr = []  # num injected @ beamoff
	num_shin_arr = []  # num shinethru @ beamoff
	for i, run in enumerate(runs):
		go = 0
		while not go:
			bdens = SimpleSignal(run, '\\bdens2_h')
			bmvol = SimpleSignal(run, '\\bmvol')
			ne = SimpleSignal(run, '\\ne')
			shinethru = SimpleSignal(run, '\\sbshine_h')  # num/sec
			sinj = SimpleSignal(run, '\\sinj')  # num/sec
			dep = SimpleSignal(run, '\\sbdepsc_h')  # \\sinj - \\sbshine_h  #dep/sec
			if sinj.data is not None and bmvol.data is not None and shinethru.data is not None: go = 1
		n0 = ne.data[0, 0]  # num/cm^3
		ntot = [np.sum(bdens.data[s, :] * bmvol.data[s, :]) for s in np.arange(len(bmvol.dim2))]
		dt = float(dep.dim1[2] - dep.dim1[1])
		num_dep = np.sum(dep.data) * dt  # total # of beam neutrals deposited during beam blip
		num_inj = np.sum(sinj.data) * dt
		num_shi = np.sum(shinethru.data) * dt
		# coupl_arr.append(1.-num_shi/num_inj)  # 1-shinethru
		
		num_lost = num_dep - np.array(
			ntot)  # number particles lost vs time, ie those deposited minus those remaining in bdens
		# ibeamoff = max(np.where(bdens.dim2 < beamon[i] + beamduration[i] + 1.e-6)[0])
		ibeamoff = np.where(num_lost == min(num_lost))[0][0]
		pr_loss_arr.append(num_lost[ibeamoff])
		coupl_arr.append(ntot[ibeamoff])  # how many remain in bdens at beamoff
		num_inj_arr.append(num_inj)
		num_shin_arr.append(num_shi)
		
		print('\n{} :: x {}'.format(run, x[i]))
		print('num_inj = {:.2e}'.format(num_inj))
		print('check 1 = {:.2f}'.format(num_lost[ibeamoff] / min(num_lost)))
		# print('check 1 = {:.2f}'.format((num_shi + num_dep) / num_inj))
		# print('check 1 = {:.2f}'.format((num_inj - num_shi) / (pr_loss_arr[-1] + coupl_arr[-1])))
		print('shinethru = {:.2f}%'.format(num_shi / num_inj * 100.))
		print('n0 = {:.2e}(cm^-3)'.format(n0))
		print('dep@beam_off = {:.2e}'.format(coupl_arr[i]))
	coupl_arr, pr_loss_arr, num_shin_arr, num_inj_arr = np.array(coupl_arr), np.array(pr_loss_arr), np.array(
		num_shin_arr), np.array(num_inj_arr)
	coupl_frac = coupl_arr / num_inj_arr
	prl_frac = pr_loss_arr / num_inj_arr
	shi_frac = num_shin_arr / num_inj_arr
	ion_frac = 1. - shi_frac
	
	if overplot:
		frac_coupled_of_deposited = coupl_arr / (num_inj_arr - num_shin_arr)
		frac_promptl_of_deposited = pr_loss_arr / (num_inj_arr - num_shin_arr)
		print('\nruns :: {}'.format(runs))
		print('check=1:: {}'.format(coupl_frac / (frac_coupled_of_deposited * ion_frac)))
		print('ionized fraction: {}'.format(ion_frac))
		print('coupled frac of depos: {}'.format(frac_coupled_of_deposited))
		print('coupled frac of injec: {}'.format(coupl_frac))
		if plot:
			ax.plot(x, frac_coupled_of_deposited, label='coupled {}'.format(label))
			ax.plot(x, frac_promptl_of_deposited, '--', label='prompt loss {}'.format(label))
	else:
		if plot:
			ax.fill_between(x, np.zeros_like(coupl_arr), coupl_frac, label='coupled {}'.format(label), alpha=0.5)
			ax.fill_between(x, coupl_frac, coupl_frac + prl_frac, label='prompt loss {}'.format(label), alpha=0.5)
			ax.fill_between(x, coupl_frac + prl_frac, coupl_frac + prl_frac + shi_frac,
			                label='shinethru {}'.format(label),
			                alpha=0.5)
			ax.plot(x, coupl_frac, 'k+-')
			ax.plot(x, coupl_frac + prl_frac, 'k+-')


class SimpleSignal:
	"""A signal class to provide an organizational scheme for storing
	timeseries data retrieved from an MDSplus signal node. To collect
	the data from an MDSplus signal node, provide the shot number,
	nodepath, and tree.
	"""
	
	def __init__(self, shot, nodepath, tree='transp_ltx'):
		"""Instatiating a signal will make an attempt to collect data
		for a particular node on the MDSplus server on
		dave.physics.wisc.edu.  Specify an integer for the shot number
		(e.g. 1161127010) and a string that describes the node path
		(e.g. '\ip'). The default tree is 'mst', but you can choose
		others available on dave.physics.wisc.edu."""
		if tree == 'mst':
			conn_name = 'dave.physics.wisc.edu'
			conn_name2 = None
			goodshot_node_test = '\ip'
		elif tree == 'ltx_b':
			conn_name = 'lithos'
			conn_name2 = None
			goodshot_node_test = '.diagnostics.magnetic.ip_rog'
		elif tree == 'transp_ltx':
			conn_name = 'transpgrid'
			conn_name2 = 'transpgrid2.pppl.gov'
			goodshot_node_test = '\\dvol'
		
		# Initialize the fields of the simpleSignal class to be the input of
		# the constructor.
		self.shot = shot
		self.nodepath = nodepath
		self.tree = tree
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
		
		# Get name of signal, last name of the nodepath
		self.name = self.nodepath.split(':')[-1].split('.')[-1]
		
		# Look at standard signal as indicator of whether or not data for this shot exists at all
		try:
			conn = MDSplus.Connection(conn_name)
			conn.openTree(tree, shot)
			testfordata = conn.get(goodshot_node_test)
			self.goodshot = True
			# Tree could have more data associated with a Signal, attempt
			# to get this information, but deal with the possibility that
			# this data may not exist for a given signal.
			try:
				node = conn.get(nodepath)
				self.data = node.data()
				self.dataunits = ' '.join(conn.get('units_of({})'.format(nodepath)).split())
				self.ndim = self.data.ndim
				self.dim1 = conn.get('dim_of({}, 0)'.format(nodepath))
				self.dim1units = ' '.join(conn.get('units_of(dim_of({}, 0))'.format(nodepath)).split())
				if self.ndim > 1:  # 2 or 3 dimensional data, get data/units for 2nd dimension
					self.dim2 = conn.get('dim_of({}, 1)'.format(nodepath))
					self.dim2units = ' '.join(conn.get('units_of(dim_of({}, 1))'.format(nodepath)).split())
				if self.ndim > 2:  # 3 dimensional data
					self.dim3 = conn.get('dim_of({}, 2)'.format(nodepath))
					self.dim3units = ' '.join(conn.get('units_of(dim_of({}, 2))'.format(nodepath)).split())
			except:
				print("{0} is not available for shot {1}".format(self.name, self.shot))
				return
		except:
			if conn_name2 is not None:
				try:
					conn = MDSplus.Connection(conn_name2)
					conn.openTree(tree, shot)
					testfordata = conn.get(goodshot_node_test)
					self.goodshot = True
					# Tree could have more data associated with a Signal, attempt
					# to get this information, but deal with the possibility that
					# this data may not exist for a given signal.
					try:
						node = conn.get(nodepath)
						self.data = node.data()
						self.dataunits = ' '.join(conn.get('units_of({})'.format(nodepath)).split())
						self.ndim = self.data.ndim
						self.dim1 = conn.get('dim_of({}, 0)'.format(nodepath))
						self.dim1units = ' '.join(conn.get('units_of(dim_of({}, 0))'.format(nodepath)).split())
						if self.ndim > 1:  # 2 or 3 dimensional data, get data/units for 2nd dimension
							self.dim2 = conn.get('dim_of({}, 1)'.format(nodepath))
							self.dim2units = ' '.join(conn.get('units_of(dim_of({}, 1))'.format(nodepath)).split())
						if self.ndim > 2:  # 3 dimensional data
							self.dim3 = conn.get('dim_of({}, 2)'.format(nodepath))
							self.dim3units = ' '.join(conn.get('units_of(dim_of({}, 2))'.format(nodepath)).split())
					except:
						print("{0} is not available for shot {1}".format(self.name, self.shot))
						return
				except:
					self.goodshot = False
					return
			else:
				self.goodshot = False
				return
		
		# Let's see if there's any more meta data like units and error bars.
		
		try:
			self.label = ' '.join(conn.get('{}.label'.format(nodepath)).split())
		except:
			pass
		if self.label is None:
			self.label = self.name
		
		# Close connection to prevent hanging processes.
		conn.closeAllTrees()
	
	def plot(self, title=None, contour=False, label=None, ax=None):
		"""Plot the signal vs. time using as much of the stored
		 information as is available to annotate the figure. Set the
		 color of the plot line using, e.g. color='b' for blue."""
		
		if ax is None:
			fig = plt.figure()
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


def smooth(x, window_len=11, window='hanning'):
	"""smooth the data using a window with requested size.

	This method is based on the convolution of a scaled window with the signal.
	The signal is prepared by introducing reflected copies of the signal
	(with the window size) in both ends so that transient parts are minimized
	in the begining and end part of the output signal.

	input:
		x: the input signal
		window_len: the dimension of the smoothing window; should be an odd integer
		window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
			flat window will produce a moving average smoothing.

	output:
		the smoothed signal

	example:

	t=linspace(-2,2,0.1)
	x=sin(t)+randn(len(t))*0.1
	y=smooth(x)

	see also:

	numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
	scipy.signal.lfilter

	TODO: the window parameter could be the window itself if an array instead of a string
	NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
	"""
	
	if x.ndim != 1:
		raise ValueError("smooth only accepts 1 dimension arrays.")
	
	if x.size < window_len:
		raise ValueError("Input vector needs to be bigger than window size.")
	
	if window_len < 3:
		return x
	
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
	
	s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
	# print(len(s))
	if window == 'flat':  # moving average
		w = np.ones(window_len, 'd')
	else:
		w = eval('np.' + window + '(window_len)')
	
	y = np.convolve(w / w.sum(), s, mode='valid')
	return y


if __name__ == '__main__':
	read_nenite('//samba/wcapecch/datasets/LTX_100981_468-1_5.nenite')

# eqdsk = 'Z:/datasets/LTX_100981_468-2_0.eqdsk'
# read_eqdsk(eqdsk)
# color_cycle_demo()
# ltx_shot = 100981
# transp_runid = 1115390512

# mst_shot = 1160926022
# mst_sig = SimpleSignal(mst_shot, '\mraw_misc::nmf', tree='mst')
# mst_sig.plot()

# get_current_nbi_shot()
# plt.subplot(311)
# ltx_sig = SimpleSignal(ltx_shot, '.oper_diags.ltx_nbi.source_diags.i_arc', tree='ltx_b')
# plt.subplot(312)
# ltx_sig.plot()
# q2d = SimpleSignal(transp_runid, '\\q', tree='transp_ltx')
# plt.subplot(313)
# q2d.plot()
# plt.show()
