try:
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt
except:
	pass
import numpy as np
import MDSplus
import re
import datetime
from bills_LTX_MDSplus_toolbox import *
import glob
import os


def is_good_shot(shot, sigs, treename='ltx_b'):
	good_shot = True
	tree = get_tree_conn(shot, treename=treename)
	for sig in sigs:
		try:
			(_, _) = get_data(tree, sig)
		except:
			good_shot = False
	return good_shot


def is_nbi_shot(shot, tree=None):
	if tree is None:
		tree = get_tree_conn(shot)
	if shot < 200000:
		i_arc_node = '.oper_diags.ltx_nbi.source_diags.i_arc'
	else:
		i_arc_node = '.ltx_nbi.source_diags.i_arc'
	times = None
	try:
		(times, dat) = get_data(tree, i_arc_node)
		nbi_shot = True
	except MDSplus.mdsExceptions.TreeNODATA:
		nbi_shot = False
		print('shot #{} is NOT a beam shot'.format(shot))
	return nbi_shot


def closest(array, value):
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	return array[idx], idx


def get_current_ltx_shot():
	# latest shot in 'ltx_b' tree
	shotdone = open('Y:\\maindata\\shotdone', 'r')
	shot = int(shotdone.read())
	
	return shot


def get_shot_timestamp(shot):
	if shot > 200000:
		treename = 'ltx_nbi'
	else:
		treename = 'ltx_b'
	tree = get_tree_conn(shot, treename=treename)
	
	return get_data(tree, '.metadata.timestamp')


def get_shots_since(date='09/01/21', nbi=True):
	since = datetime.datetime.strptime(date, '%m/%d/%y')
	
	lastshot = get_current_ltx_shot()
	shot_arr = []
	while 1:
		# print(lastshot)
		try:
			t = get_tree_conn(lastshot, treename='ltx_b')
			ts = get_data(t, '.metadata.timestamp')
			print(f'shot:{lastshot}, time:{ts}')
			dt = datetime.datetime.strptime(ts, '%m/%d/%y %I:%M:%S %p')
			if dt < since:
				break
			
			if nbi:
				if is_nbi_shot(lastshot, t):
					shot_arr.append(lastshot)
			else:
				shot_arr.append(lastshot)
		except MDSplus.mdsExceptions.TreeFILE_NOT_FOUND:
			print(f'no tree date for shot {lastshot}')
			pass
		lastshot -= 1
	return shot_arr


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
	# THIS VERSION assumes each row of psixy start on newline- not true at least for TRANSP output GEQDSK files- use read_eqdsk3
	f = open(eqdsk, 'r')
	form = r'\S\d+.\d+E[+-]\d+'
	_, nxeqd, nyeqd = [int(dig) for dig in re.findall(r'\d+', f.readline())][
	                  -3:]  # transp version has more numbers on leadin
	xdimeqd, ydimeqd, r0, redeqd, ymideqd = [float(n) for n in re.findall(form, f.readline())]
	xma, yma, psimag, psilim, beqd = [float(n) for n in re.findall(form, f.readline())]
	toteqd, psimx1, psimx2, xax1, xax2 = [float(n) for n in re.findall(form, f.readline())]
	zax1, zax2, psisep, xsep, ysep = [float(n) for n in re.findall(form, f.readline())]
	nrows = np.ceil(nxeqd / 5.)  # num rows to scan fpol, etc
	nrows_psi = np.ceil(nxeqd * nyeqd / 5.)  # num rows to scan for psixy
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
	
	rvec_hd, zvec_hd = np.linspace(x[0], x[-1], num=1000), np.linspace(y[0], y[-1], num=1000)
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
	          'rves': rves, 'zves': zves, 'x': x, 'y': y, 'psi': psiout, 'psixy': psixy, 'real_psilim': real_psilim,
	          'q_xy': q_xy,
	          'bphi_xy': bphi_xy, 'bp_xy': bp_xy, 'x_xy': x_xy, 'y_xy': y_xy, 'br_xy': br_xy, 'bz_xy': bz_xy,
	          'bp_mp': bp_mp, 'bphi_mp': bphi_mp, 'bave_mp': bave_mp, 'bave_xy': bave_xy, 'br_mp': br_mp,
	          'bz_mp': bz_mp, 'ip_is_cw': ip_is_cw, 'bt_is_cw': bt_is_cw}
	
	# fpol_xy: fpol_xy, rminor_xy: rminor_xy, norm_psixy: norm_psixy, plasma: new_plasma, dpsi: dpsi,
	# psiaxis: psi_axis}
	return output


def read_eqdsk3(eqdsk, plot=False):
	# THIS VERSION assumes psixy data is flattened so new rows of psixy do not go onto newline in data (as in TRANSP GEQDSK output).
	# For data that does put new psixy rows onto newlines, use read_eqdsk().
	f = open(eqdsk, 'r')
	form = r'\S\d+.\d+E[+-]\d+'
	_, nxeqd, nyeqd = [int(dig) for dig in re.findall(r'\d+', f.readline())][
	                  -3:]  # transp version has more numbers on leadin
	xdimeqd, ydimeqd, r0, redeqd, ymideqd = [float(n) for n in re.findall(form, f.readline())]
	xma, yma, psimag, psilim, beqd = [float(n) for n in re.findall(form, f.readline())]
	toteqd, psimx1, psimx2, xax1, xax2 = [float(n) for n in re.findall(form, f.readline())]
	zax1, zax2, psisep, xsep, ysep = [float(n) for n in re.findall(form, f.readline())]
	nrows = np.ceil(nxeqd / 5.)  # num rows to scan fpol, etc
	nrows_psi = np.ceil(nxeqd * nyeqd / 5.)  # num rows to scan for psixy
	fpol, pres, ffpeqd, ppeqd, psixy, q = [], [], [], [], [], []
	for i in np.arange(nrows) + 1:
		fpol.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		pres.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		ffpeqd.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows) + 1:
		ppeqd.extend([float(n) for n in re.findall(form, f.readline())])
	for i in np.arange(nrows_psi) + 1:
		psixy.extend([float(n) for n in re.findall(form, f.readline())])
	psixy = np.array(psixy).reshape((nyeqd, nxeqd))
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
	
	x = np.linspace(redeqd, redeqd + xdimeqd, endpoint=True, num=nxeqd)
	y = np.linspace(ymideqd - ydimeqd / 2., ymideqd + ydimeqd / 2., endpoint=True, num=nyeqd)
	[x_xy, y_xy] = np.meshgrid(x, y)
	
	psiout = np.linspace(psimag, psilim, nxeqd)
	real_psilim = (psilim - .01 * psimag) / .99
	psirzvec = psixy.reshape(nxeqd * nyeqd)
	q_xy = np.interp(psirzvec, psiout, q)
	BtR = np.interp(psirzvec, psiout, fpol)
	BtR = BtR.reshape([nyeqd, nxeqd])
	# bphi_xy = BtR / np.transpose(x_xy)  # WJC commented out 6/23/21- why the transpose?? looks wrong
	bphi_xy = BtR / x_xy
	bphi_mp = fpol / x
	
	rvec_hd, zvec_hd = np.linspace(x[0], x[-1], num=1000), np.linspace(y[0], y[-1], num=1000)
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
		lr, lz, _, _ = ltx_limiter()
		fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey='row')
		ax1.contourf(x_xy, y_xy, br_xy)
		ax1.plot(rlimit, zlimit, 'k-')
		ax2.contourf(x_xy, y_xy, bz_xy)
		ax3.contourf(x_xy, y_xy, bphi_xy)
		for ax in [ax1, ax2, ax3]:
			ax.plot(lr, lz, 'k-')
		ax1.set_title('Br')
		ax2.set_title('Bz')
		ax3.set_title('Bphi')
		# plt.plot(x, br_xy[nz0, :], label='br_xy')
		# plt.plot(x, bz_xy[nz0, :], label='bz_xy')
		# plt.plot(x, bphi_mp, label='bt')
		# plt.legend()
		fig2, (ax4, ax5) = plt.subplots(ncols=2)
		ax4.contour(x_xy, y_xy, psixy, levels=np.linspace(psimag, psilim, num=15), linestyles='--', colors='k')
		ax4.plot(lr, lz, 'k-')
		ax4.plot(rlimit, zlimit, 'r-')
		ax4.axis('equal')
		ax4.set_ylim((-1.1 * max(lz), 1.1 * max(lz)))
		ax4.axis('off')
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
	          'rves': rves, 'zves': zves, 'x': x, 'y': y, 'psi': psiout, 'psixy': psixy, 'real_psilim': real_psilim,
	          'q_xy': q_xy,
	          'bphi_xy': bphi_xy, 'bp_xy': bp_xy, 'x_xy': x_xy, 'y_xy': y_xy, 'br_xy': br_xy, 'bz_xy': bz_xy,
	          'bp_mp': bp_mp, 'bphi_mp': bphi_mp, 'bave_mp': bave_mp, 'bave_xy': bave_xy, 'br_mp': br_mp,
	          'bz_mp': bz_mp, 'ip_is_cw': ip_is_cw, 'bt_is_cw': bt_is_cw}
	
	# fpol_xy: fpol_xy, rminor_xy: rminor_xy, norm_psixy: norm_psixy, plasma: new_plasma, dpsi: dpsi,
	# psiaxis: psi_axis}
	return output


def read_eqdsk2(filename, plot=False):
	def read_1d(fid, j, n):
		output = np.zeros((n,))
		for i in range(n):
			if j == 0:
				line = fid.readline()
			output[i] = line[j:j + 16]
			j += 16
			if j == 16 * 5:
				j = 0
		return output, j
	
	def read_2d(fid, j, n, m):
		output = np.zeros((m, n))
		for k in range(n):
			for i in range(m):
				if j == 0:
					line = fid.readline()
				output[i, k] = line[j:j + 16]
				j += 16
				if j == 16 * 5:
					j = 0
		return output, j
	
	# Read-in data
	eqdsk_obj = {}
	with open(filename, 'r') as fid:
		# Get sizes
		line = fid.readline()
		split_line = line.split()
		eqdsk_obj['nz'] = int(split_line[-1])
		eqdsk_obj['nr'] = int(split_line[-2])
		# Read header content
		line_keys = [['rdim', 'zdim', 'raxis', 'rleft', 'zmid'],
		             ['raxis', 'zaxis', 'psimax', 'psimin', 'bcentr'],
		             ['itor', 'skip', 'skip', 'skip', 'skip'],
		             ['skip', 'skip', 'skip', 'skip', 'skip']]
		for i in range(4):
			line = fid.readline()
			for j in range(5):
				if line_keys[i][j] == 'skip':
					continue
				line_seg = line[j * 16:(j + 1) * 16]
				eqdsk_obj[line_keys[i][j]] = float(line_seg)
		# Read flux profiles
		j = 0
		keys = ['fpol', 'pres', 'ffprime', 'pprime']
		for key in keys:
			eqdsk_obj[key], j = read_1d(fid, j, eqdsk_obj['nr'])
		# Read PSI grid
		eqdsk_obj['psirz'], j = read_2d(fid, j, eqdsk_obj['nz'],
		                                eqdsk_obj['nr'])
		# Read q-profile
		eqdsk_obj['qpsi'], j = read_1d(fid, j, eqdsk_obj['nr'])
		# Skip line (data already present)
		line = fid.readline()
		# Read outer flux surface
		eqdsk_obj['rzout'], j = read_2d(fid, j, eqdsk_obj['nr'], 2)
		# Read limiting corners
		eqdsk_obj['rzlim'], j = read_2d(fid, j, 5, 2)
	return eqdsk_obj


def ltx_limiter(plot=False):
	'''
	myFile = fopen('ltx_limiter.txt','r');
	myNum = fscanf(myFile,'%d',1);
	rlimiter = fscanf(myFile,'%f',myNum);
	zlimiter = fscanf(myFile,'%f',myNum);
	nlim = numel(rlimiter);
	fclose(myFile);
	'''
	lim_fn = 'Z:/PycharmProjects/LTXb-py/ltx_limiter.txt'
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
	
	# adding in pts along gaps on LF and HFS
	rbot, zbot = rlimiter[:103], zlimiter[:103]
	rtop, ztop = rlimiter[103:], zlimiter[103:]
	nh, nl = int((ztop[-1] - zbot[0]) * 100.), int((ztop[0] - zbot[-1]) * 100.)
	rhfs, zhfs = np.ones(nh) * min(rlimiter), np.linspace(ztop[-1], zbot[0], endpoint=True, num=nh + 2)
	zhfs = zhfs[1:-1]
	rpos, zpos = rhfs[np.where(zhfs > 0)], zhfs[np.where(zhfs > 0)]
	rneg, zneg = rhfs[np.where(zhfs < 0)], zhfs[np.where(zhfs < 0)]
	rlfs, zlfs = np.ones(nl) * max(rlimiter), np.linspace(zbot[-1], ztop[0], endpoint=True, num=nl + 2)
	zlfs = zlfs[1:-1]
	rlimiter, zlimiter = np.array([]), np.array([])
	for arr in [[min(rbot)], rneg, rbot, rlfs, rtop, rpos, [min(rbot)]]:
		rlimiter = np.append(rlimiter, arr)
	for arr in [[0.], zneg, zbot, zlfs, ztop, zpos, [0.]]:
		zlimiter = np.append(zlimiter, arr)

	rmag = 0.4  # doesn't need to be exact, just using interior pt. to find pts outside limiter
	rminor_lim = np.sqrt((rlimiter - rmag) ** 2 + zlimiter ** 2)
	theta_lim = np.arctan2(zlimiter, (rlimiter - rmag))
	theta_lim[0] -= 2. * np.pi
	# theta_lim[-1] += 2. * np.pi
	
	if plot:
		plt.plot(rlimiter, zlimiter)
		plt.axis('equal')
		plt.show()
	
	return rlimiter, zlimiter, rminor_lim, theta_lim


def ltx_above():
	rl, _, _, _ = ltx_limiter()
	th = np.linspace(-np.pi, np.pi, endpoint=True, num=250)
	xin, xout, yin, yout = min(rl) * np.cos(th), max(rl) * np.cos(th), min(rl) * np.sin(th), max(rl) * np.sin(th)
	return xin, xout, yin, yout


def make_patch_spines_invisible(ax):
	ax.set_frame_on(True)
	ax.patch.set_visible(False)
	for sp in ax.spines.values():
		sp.set_visible(False)


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


def calculate_perv_pwr(shot):
	if shot > 200000:
		treename = 'ltx_nbi'
		prefix = ''
		nbi_only = True
	else:
		treename = 'ltx_b'
		prefix = '.oper_diags.ltx_nbi'
		nbi_only = False
	try:
		tree = get_tree_conn(shot, treename=treename)
		(ti, ib) = get_data(tree, f'{prefix}.source_diags.i_hvps')
		(tv, vb) = get_data(tree, f'{prefix}.source_diags.v_hvps')
		ib = np.interp(tv, ti, ib)  # get ib onto vb axis
		t_beamon = np.where(vb > 5000)  # only look at where beam is above 5kV
		pad = 0.25e-3  # remove this amt from beginning/end of beam
		t_window = np.where((tv >= tv[t_beamon[0][0]] + pad) & (tv <= tv[t_beamon[0][-1]] - pad))
		perv, pwr = np.ones_like(vb), np.ones_like(vb)
		perv[:], pwr[:] = np.nan, np.nan
		perv[t_window] = ib[t_window] / vb[t_window] ** 1.5
		pwr[t_window] = ib[t_window] * vb[t_window]  # [W]
		av_perv, av_pwr = np.nanmean(perv), np.nanmean(pwr)
		dperv, dpwr = (np.nanmax(perv) - np.nanmin(perv)) / 2., (np.nanmax(pwr) - np.nanmin(pwr)) / 2.
		return av_perv, dperv, av_pwr, dpwr
	except mds.mdsExceptions.TreeNODATA:
		print(f'trouble fetching MDSPlus data for shot {shot}')
		return np.nan, np.nan, np.nan, np.nan


def avg_perv(shot, Ej=False, ib_return=False, twin=None):
	if shot > 200000:
		treename = 'ltx_nbi'
		prefix = ''
	else:
		treename = 'ltx_b'
		prefix = '.oper_diags.ltx_nbi'
	tree = get_tree_conn(shot, treename=treename)
	(ti, ib) = get_data(tree, f'{prefix}.source_diags.i_hvps')
	(tv, vb) = get_data(tree, f'{prefix}.source_diags.v_hvps')
	ib = np.interp(tv, ti, ib)  # get ib onto vb axis
	# t_beamon = np.where(vb > 5000)  # only look at where beam is above 5kV
	if twin is None:
		# t_beamon = np.where(
		# 	(.46 < tv) & (tv < .463))  # look between 460-463ms for these shots (ignore rough Ibeam startup)
		ibeamon = np.where(vb > 5000)[0]
		dt = 0.25e-3
		(ton, toff) = (tv[ibeamon[0]] + dt, tv[ibeamon[-1]] - dt)
		t_beamon = np.where((ton < tv) & (tv < toff))
	else:
		t_beamon = np.where((twin[0] < tv) & (tv < twin[1]))
	pb = ib * vb  # beam power [W]
	tot_joules = np.sum([pb[i] * (tv[i] - tv[i - 1]) for i in np.arange(1, len(tv))])
	if len(t_beamon[0]) < 10:
		av_perv = np.nan
		av_ib = np.nan
	else:
		perv = np.ones_like(vb)
		perv[:] = np.nan
		perv[t_beamon] = ib[t_beamon] / vb[t_beamon] ** 1.5
		av_perv = np.nanmean(perv)
		av_ib = np.nanmean(ib[t_beamon])
	if Ej:
		if ib_return:
			return av_perv, tot_joules, av_ib
		else:
			return av_perv, tot_joules
	if ib_return:
		return av_perv, av_ib
	else:
		return av_perv


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


def multiple_legends_example():
	d1 = np.linspace(0, 20)
	d2, d3 = np.sin(d1), np.cos(d1)
	fig, ax = plt.subplots()
	p1, = ax.plot(d1, d2, label='sin')
	p2, = ax.plot(d1, d3, label='cos')
	p3, = ax.plot(np.nan, np.nan, label=f'dude')
	p4, = ax.plot(np.nan, np.nan, label=f'bro')
	leg1 = ax.legend(handles=[p1, p2], loc='upper left', title='whoa')
	ax.legend(handles=[p3, p4], loc='center left')
	ax.add_artist(leg1)
	
	plt.show()


def read_transp_cdf(run, local_dir='//samba/wcapecch/transp_rawdata/'):
	# NOTE this is poorly outfitted- requires looking in /12/ directory which may not be the case
	import os
	from scipy.io import netcdf
	print("Looking for CDF data for run {}... be patient".format(run))
	# convert numeric into alpha-numeric
	cdf_fn = str(run)[:6] + chr(ord('@') + int(str(run)[6:8])) + str(run)[8:] + '.CDF'
	if os.path.isdir('//p/transparch/result/LTX/12/'):
		dir = '//p/transparch/result/LTX/12/'
	else:
		print(f'Running locally... searching {local_dir}')
	fp = '{}{}'.format(local_dir, cdf_fn)
	if os.path.isfile(fp):
		cdf = netcdf.netcdf_file(fp, 'r', mmap=False).variables
		return cdf
	else:
		print('TRANSP cdf {} NOT FOUND in {}'.format(cdf_fn, dir))
		return None


def smooth(x, window_len=11, window='hanning', mode='valid'):
	"""smooth the data using a window with requested size.

	This method is based on the convolution of a scaled window with the signal.
	The signal is prepared by introducing reflected copies of the signal
	(with the window size) in both ends so that transient parts are minimized
	in the beginning and end part of the output signal.

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
	
	if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError("Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
	
	if mode not in ['valid', 'same']:
		raise ValueError("Mode is one of 'valid', 'same'")
	
	if mode is 'valid':
		s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
	else:
		s = x
	# print(len(s))
	if window == 'flat':  # moving average
		w = np.ones(window_len, 'd')
	else:
		w = eval('np.' + window + '(window_len)')
	
	if mode is 'valid':
		y = np.convolve(w / w.sum(), s, mode='valid')
		y = y[int(window_len / 2):-int(window_len / 2)]  # undo addition in line 808
	else:
		y = np.convolve(w / w.sum(), s, mode='same')
	return y


if __name__ == '__main__':
	eq = 'Z:/transp/t106536/106536R02_05.eqdsk'
	read_eqdsk3(eq, plot=True)
	# multiple_legends_example()
	a = 1
# ll = get_shots_since(date='09/14/21')
# a = 1
# ltx_limiter()
# read_nenite('//samba/wcapecch/datasets/LTX_100981_468-1_5.nenite')

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
