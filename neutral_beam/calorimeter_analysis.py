import matplotlib.gridspec as gridspec
from toolbox.helpful_stuff import avg_perv, calculate_perv_pwr, get_shot_timestamp
import lvm_read
from bills_LTX_MDSplus_toolbox import *
import os
from scipy.optimize import minimize

plot_shot_signals = 0
if plot_shot_signals:
	shot = 104253
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	fig, ax = plt.subplots()
	tree = get_tree_conn(shot, treename='ltx_b')
	prefix = '.oper_diags.ltx_nbi'
	(ti, ib) = get_data(tree, f'{prefix}.source_diags.i_hvps')
	(tv, vb) = get_data(tree, f'{prefix}.source_diags.v_hvps')
	ax.plot(tv, vb / 1000., c=clrs[0])
	axr = ax.twinx()
	axr.plot(ti, ib, c=clrs[1])
	fs = 12
	ax.set_xlabel('time (s)', fontsize=fs)
	ax.set_ylabel('$V_{nbi}$ (kV)', c=clrs[0], fontsize=fs)
	axr.set_ylabel('$I_{nbi}$ (A)', c=clrs[1], fontsize=fs)
	plt.show()


def cal_dtemp(shot, dbug=False, more=False):
	if more:
		nodata = np.nan, np.nan, [np.nan, np.nan, np.nan, np.nan, np.nan], [np.nan, np.nan]
	else:
		nodata = np.nan, np.nan
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
		ts = get_data(tree, '.metadata.timestamp')
	except mds.mdsExceptions.TreeNODATA:
		print(f'trouble fetching MDSPlus data for shot {shot}')
		return nodata
	if len(ts.split('/')[0]) < 2:
		ts = f'0{ts}'
	direc = 'Y:/thermocouple/Calorimeter/'
	if nbi_only:
		dt = datetime.datetime.strptime(ts, '%m/%d/%Y %I:%M:%S %p')
		lvm_files = os.listdir(direc)
		lvm_files = [lvmf for lvmf in lvm_files if lvmf.endswith('.lvm') and len(lvmf) == 18]  # lvms with date format
		lvm_times = [datetime.datetime.strptime(lvmf.split('.')[0], '%m%d%Y%H%M%S') for lvmf in lvm_files]
		sync_file = [lvm_files[i] for i in range(len(lvm_files)) if
		             dt < lvm_times[i] and dt + datetime.timedelta(seconds=60) > lvm_times[i]]
	else:
		if os.path.exists(f'{direc}{shot}.lvm'):
			sync_file = [f'{shot}.lvm']
		else:
			print(f'---- lvm file not found for shot {shot}')
			return nodata
	if len(sync_file) != 1:
		print(f'---- problem syncing shot {shot}, timestamp {ts} with lvm file')
		return nodata
	else:
		print(f'shot {shot} synced with {sync_file[0]}')
		ib = np.interp(tv, ti, ib)  # get ib onto vb axis
		t_beamon = np.where(vb > 5000)  # only look at where beam is above 5kV
		pad = 0.25e-3  # remove this amt from beginning/end of beam
		t_window = np.where((tv >= tv[t_beamon[0][0]] + pad) & (tv <= tv[t_beamon[0][-1]] - pad))
		ipfit = np.polyfit(tv[t_window], ib[t_window], 1)
		stdev = np.sqrt(np.sum((ib[t_window] - (tv[t_window] * ipfit[0] + ipfit[1])) ** 2 / len(tv[t_window])))
		pb = ib * vb  # beam power [W]
		tot_joules = np.sum([pb[i] * (tv[i] - tv[i - 1]) for i in np.arange(1, len(tv))])
		
		lvm_file = sync_file[0]
		lvm = lvm_read.read(f'{direc}{lvm_file}', dump_file=False, read_from_pickle=False)
		time, rtd_nam, c0, c1, c2, c3, c4 = get_cal_sigs(lvm[0])
		it = np.where(time < 300)
		time, c0, c1, c2, c3, c4 = time[it], c0[it], c1[it], c2[it], c3[it], c4[it]
		ipeak = np.where(c0 == max(c0))[0][0]
		baseline = mean((c0[:ipeak - 1] + c1[:ipeak - 1] + c2[:ipeak - 1] + c3[:ipeak - 1] + c4[:ipeak - 1]) / 5.)
		islow_decay = np.where((time > time[ipeak] + 25) & (time < time[ipeak] + 100))[0]
		if len(islow_decay) < 50:
			print(f'--insufficient data to analyze for shot {shot} in file {lvm_file}')
			return nodata
		else:
			dt0, dt1, dt2, dt3, dt4 = c0[ipeak] - mean(c0[0:ipeak - 5]), c1[ipeak + 5] - mean(c1[0:ipeak - 5]), c2[
				ipeak + 5] - mean(c2[0:ipeak - 5]), c3[ipeak + 5] - mean(c3[0:ipeak - 5]), c4[ipeak + 5] - mean(
				c4[0:ipeak - 5])
			line2fit = (c0[islow_decay] + c1[islow_decay] + c2[islow_decay] + c3[islow_decay] + c4[islow_decay]) / 5.
			time2fit = time[islow_decay]
			fit = np.polyfit(time2fit, line2fit, 1)
			temp_interp = np.interp(time[ipeak], time, time * fit[0] + fit[1])
			meas = temp_interp - baseline
			pred = tot_joules / c_ps / m_copp  # predicted dT (degC) from injected beam energy
			
			if dbug:
				plt.suptitle(f'{shot}')
				for i, c in enumerate([c0, c1, c2, c3, c4]):
					plt.plot(time, c, label=rtd_nam[i])
				# plt.axvline(time[islow_decay[0]])
				# plt.axvline(time[islow_decay[-1]])
				plt.legend()
				plt.xlim((0, 60))
				# plt.axhline(baseline)
				# plt.plot(time, time * fit[0] + fit[1])
				plt.xlabel('time (s)')
				plt.ylabel('temp ($\degree$C)')
				plt.show()
			if more:
				return meas, pred, [dt0, dt1, dt2, dt3, dt4], [ipfit[0], stdev]
			else:
				return meas, pred


def get_cal_sigs(lvm):
	time = lvm['data'][:, 0]
	if 'Temperature_0' in lvm['Channel names']:  # data taken using fast calorimeter-only LabVIEW
		rtd_nam = lvm['Channel names'][4:9]  # tc14-18
		c0, c1, c2 = lvm['data'][:, 4], lvm['data'][:, 5], lvm['data'][:, 6]
		c3, c4 = lvm['data'][:, 7], lvm['data'][:, 8]
	else:  # data taken with regular LabVIEW, records all RTDs along with calorimeter thermocouples
		rtd_nam = lvm['Channel names'][15:20]  # tc14-18
		c0, c1, c2 = lvm['data'][:, 15], lvm['data'][:, 16], lvm['data'][:, 17]
		c3, c4 = lvm['data'][:, 18], lvm['data'][:, 19]
	return time, rtd_nam, c0, c1, c2, c3, c4


def cal_temp_sigs(shot, calonly=False):
	nodata = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
	if calonly:
		nbi_only = False
	else:
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
			ts = get_data(tree, '.metadata.timestamp')
		except mds.mdsExceptions.TreeNODATA:
			print(f'trouble fetching MDSPlus data for shot {shot}')
			return nodata
		if len(ts.split('/')[0]) < 2:
			ts = f'0{ts}'
	direc = 'Y:/thermocouple/Calorimeter/'
	if nbi_only:
		dt = datetime.datetime.strptime(ts, '%m/%d/%Y %I:%M:%S %p')
		lvm_files = os.listdir(direc)
		lvm_files = [lvmf for lvmf in lvm_files if lvmf.endswith('.lvm') and len(lvmf) == 18]  # lvms with date format
		lvm_times = [datetime.datetime.strptime(lvmf.split('.')[0], '%m%d%Y%H%M%S') for lvmf in lvm_files]
		sync_file = [lvm_files[i] for i in range(len(lvm_files)) if
		             dt < lvm_times[i] and dt + datetime.timedelta(seconds=60) > lvm_times[i]]
	else:
		if os.path.exists(f'{direc}{shot}.lvm'):
			sync_file = [f'{shot}.lvm']
		else:
			print(f'---- lvm file not found for shot {shot}')
			return nodata
	if len(sync_file) != 1:
		print(f'---- problem syncing shot {shot}, timestamp {ts} with lvm file')
		return nodata
	else:
		print(f'shot {shot} synced with {sync_file[0]}')
		lvm_file = sync_file[0]
		lvm = lvm_read.read(f'{direc}{lvm_file}', dump_file=False, read_from_pickle=False)
		time, rtd_nam, c0, c1, c2, c3, c4 = get_cal_sigs(lvm[0])
		it = np.where(time < 300)
		time, c0, c1, c2, c3, c4 = time[it], c0[it], c1[it], c2[it], c3[it], c4[it]
		ipeak = min([np.where(c0 == max(c0))[0][0], np.where(c1 == max(c1))[0][0], np.where(c2 == max(c2))[0][0],
		             np.where(c3 == max(c3))[0][0], np.where(c4 == max(c4))[0][0]])
		b0, b1, b2, b3, b4 = mean(c0[:ipeak - 5]), mean(c1[:ipeak - 5]), mean(c2[:ipeak - 5]), mean(
			c3[:ipeak - 5]), mean(c4[:ipeak - 5])
		c0, c1, c2, c3, c4 = c0 - b0, c1 - b1, c2 - b2, c3 - b3, c4 - b4
		time -= time[ipeak]
		return time, c0, c1, c2, c3, c4


c_ps = 385.  # specific heat of copper J/kg/degC
m_copp = 1.5834  # calorimeter copper mass [kg]


def big_scan():
	bad = [506588, 506597, 506602]
	shots11Oct = np.array(
		[506587, 506588, 506589, 506590, 506591, 506592, 506593, 506594, 506596, 506597, 506598, 506600, 506601, 506602,
		 506603, 506609, 506610, 506612, 506613, 506614, 506616, 506617, 506619, 506620, 506621, 506622, 506626, 506627,
		 506628, 506630, 506632, 506633, 506634, 506636, 506637, 506638, 506639, 506642, 506643, 506644, 506646, 506647,
		 506648])
	i60c, i60v, i70, i80 = np.where(shots11Oct < 506610)[0], np.where((shots11Oct >= 506610) & (shots11Oct < 506633))[
		0], \
	                       np.where((shots11Oct >= 506633) & (shots11Oct < 506640))[0], np.where(shots11Oct >= 506640)[
		                       0]
	dt_pred_60psi, dt_meas_60psi = [], []  # shot arrays for joules
	for ish, shot in enumerate(shots11Oct):
		dtm, dtp = cal_dtemp(shot, dbug=True)
		dt_meas_60psi.append(dtm)
		dt_pred_60psi.append(dtp)
	
	fig, ax1 = plt.subplots()  # , figsize=(12, 6))
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	slope60c = np.nanmean([dt_meas_60psi[i] / dt_pred_60psi[i] for i in i60c])
	slope60v = dt_meas_60psi[np.where(shots11Oct == 506632)[0][0]] / dt_pred_60psi[np.where(shots11Oct == 506632)[0][0]]
	slope70 = dt_meas_60psi[np.where(shots11Oct == 506639)[0][0]] / dt_pred_60psi[np.where(shots11Oct == 506639)[0][0]]
	slope80 = dt_meas_60psi[np.where(shots11Oct == 506648)[0][0]] / dt_pred_60psi[np.where(shots11Oct == 506648)[0][0]]
	
	ax1.plot([0, max(dt_pred_60psi)], [0, max(dt_pred_60psi) * slope60c], '--', c=clrs[0])
	ax1.plot(np.array(dt_pred_60psi)[i60c], np.array(dt_meas_60psi)[i60c], 'o', label='60psi same neutralizer fill',
	         c=clrs[0])
	ax1.annotate(f'{slope60c * 100.:.2f}%', (max(dt_pred_60psi), max(dt_pred_60psi) * slope60c), c=clrs[0])
	
	ax1.plot([0, max(dt_pred_60psi)], [0, max(dt_pred_60psi) * slope60v], '--', c=clrs[1])
	ax1.plot(np.array(dt_pred_60psi)[i60v], np.array(dt_meas_60psi)[i60v], 'o', label='60psi varying neutralizer',
	         c=clrs[1])
	ax1.annotate(f'{slope60v * 100.:.2f}%', (max(dt_pred_60psi), max(dt_pred_60psi) * slope60v), c=clrs[1])
	
	ax1.plot([0, max(dt_pred_60psi)], [0, max(dt_pred_60psi) * slope70], '--', c=clrs[2])
	ax1.plot(np.array(dt_pred_60psi)[i70], np.array(dt_meas_60psi)[i70], 'o', label='70psi varying neutralizer',
	         c=clrs[2])
	ax1.annotate(f'{slope70 * 100.:.2f}%', (max(dt_pred_60psi), max(dt_pred_60psi) * slope70), c=clrs[2])
	
	ax1.plot([0, max(dt_pred_60psi)], [0, max(dt_pred_60psi) * slope80], '--', c=clrs[3])
	ax1.plot(np.array(dt_pred_60psi)[i80], np.array(dt_meas_60psi)[i80], 'o', label='80psi varying neutralizer',
	         c=clrs[3])
	ax1.annotate(f'{slope80 * 100.:.2f}%', (max(dt_pred_60psi), max(dt_pred_60psi) * slope80), c=clrs[3])
	
	shots13Oct = 104000 + np.array([55, 56, 58, 60, 61, 62, 63, 64, 66, 67, 68, 69, 71, 74, 75, 76, 77, 78, 79])
	dtbeam = np.array([5, 5, 5, 5, 6, 6, 6, 6, 6, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2])
	dt_pred_60psi, dt_meas_60psi = [], []  # shot arrays for joules
	i6ms, i4ms, i2ms = np.where(dtbeam == 6), np.where(dtbeam == 4), np.where(dtbeam == 2)
	for ish, shot in enumerate(shots13Oct):
		dtm, dtp = cal_dtemp(shot)
		dt_meas_60psi.append(dtm)
		dt_pred_60psi.append(dtp)
	ax1.plot(np.array(dt_pred_60psi)[i6ms], np.array(dt_meas_60psi)[i6ms], 'o', label='60psi 6ms pulse', c=clrs[4])
	ax1.plot(np.array(dt_pred_60psi)[i4ms], np.array(dt_meas_60psi)[i4ms], 'o', label='60psi 4ms pulse', c=clrs[5])
	ax1.plot(np.array(dt_pred_60psi)[i2ms], np.array(dt_meas_60psi)[i2ms], 'o', label='60psi 2ms pulse', c=clrs[6])
	
	# ax1.plot(np.linspace(0, max(dt_pred)), np.linspace(0, max(dt_pred)), 'k--', alpha=.5)
	ax1.set_xlabel('predicted $\Delta T$ ($\degree$C)')
	ax1.set_ylabel('measured $\Delta T$ ($\degree$C)')
	ax1.set_title('Calorimeter')
	ax1.legend(fontsize=10)


def pulse_timing():
	shots60psi = 104000 + np.array([83, 84, 87, 91, 88, 92, 94, 95, 96, 97, 98, 99, 100, 101])
	nominal_strt_60psi = [-6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
	beam_strt_rel_to_neut_puff_60psi = [nom + 5 for nom in nominal_strt_60psi]
	
	shots100psi = 104000 + np.array([107, 108, 110, 106, 115, 117, 121, 124, 127, 131, 133])
	nominal_strt_100psi = [-6, -4, -2, 0, 2, 4, 6, 8, 12, 16, 20]
	beam_strt_rel_to_neut_puff_100psi = [nom + 5 for nom in nominal_strt_100psi]
	
	shots100psi_ext = 104000 + np.array([109, 111, 114, 116, 120, 122, 125, 128, 132, 134])
	nominal_strt_100psi_ext = [-4, -2, 0, 2, 4, 6, 8, 12, 16, 20]
	beam_strt_rel_to_neut_puff_100psi_ext = [nom + 5 for nom in nominal_strt_100psi_ext]
	
	dt_pred_60psi, dt_meas_60psi = [], []
	dt_pred_100psi, dt_meas_100psi = [], []
	dt_pred_100psi_ext, dt_meas_100psi_ext = [], []
	for ish, shot in enumerate(shots60psi):
		dtm, dtp = cal_dtemp(shot)
		dt_meas_60psi.append(dtm)
		dt_pred_60psi.append(dtp)
	for ish, shot in enumerate(shots100psi):
		dtm, dtp = cal_dtemp(shot)
		dt_meas_100psi.append(dtm)
		dt_pred_100psi.append(dtp)
	for ish, shot in enumerate(shots100psi_ext):
		dtm, dtp = cal_dtemp(shot)
		dt_meas_100psi_ext.append(dtm)
		dt_pred_100psi_ext.append(dtp)
	
	fig, ax1 = plt.subplots()  # , figsize=(12, 6))
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	ax1.plot(beam_strt_rel_to_neut_puff_60psi, np.array(dt_meas_60psi) / np.array(dt_pred_60psi), 'o-', c=clrs[0],
	         label='5ms 60psi')
	ax1.plot(beam_strt_rel_to_neut_puff_100psi, np.array(dt_meas_100psi) / np.array(dt_pred_100psi), 'o-', c=clrs[1],
	         label='5ms 100psi puff')
	ax1.plot(beam_strt_rel_to_neut_puff_100psi_ext, np.array(dt_meas_100psi_ext) / np.array(dt_pred_100psi_ext), 'o-',
	         label='5ms+ 100psi puff', c=clrs[2])
	ax1.axvline(5, ls='--', c='k', alpha=0.5)
	ax1.annotate('nominal', (5, .5))
	ax1.set_xlabel('tbeam_start - tpuff_start (ms)')
	ax1.set_ylabel('predicted/measured dT')
	ax1.set_title('Calorimeter')
	ax1.legend(fontsize=10)
	print(f'max fraction: {np.max(np.array(dt_meas_100psi_ext) / np.array(dt_pred_100psi_ext))}')


def anode_cathode_puffing():
	# data taken 10/15/21, 10/19/21
	# adjusting anode/cathode puff timings to see if grid filling is the issue
	
	def plot_some_data(duration, shots, clr, lbl, ax):
		if len(shots) > 0:
			dt_pred, dt_meas = [], []
			for ish, shot in enumerate(shots):
				dtm, dtp = cal_dtemp(shot)
				dt_meas.append(dtm)
				dt_pred.append(dtp)
			temp_fraction = np.array(dt_meas) / np.array(dt_pred)
			ax.plot(duration, temp_fraction, 'o', c=clr, label=lbl)
		print(f'max {max(temp_fraction)}')
	
	# 100psi on neutralizer, 60 on cathode, 30 on anode
	shots_anode = 104000 + np.array([177, 179, 182, 183, 184, 186, 187, 189, 192])
	anode_strt = np.array([7, 5, 3, 3, 0, 7, 7, 7, 7])
	anode_duration = np.array([2, 4, 6, 6, 9, 4, 6, 8, 0])
	shots_anode_noneut = 104000 + np.array([185, 191])
	anode_strt_noneut = np.array([0, 7])
	anode_duration_noneut = np.array([9, 8])
	
	# GVC, GVA = 60psi, GVN = 100psi
	shots_anode2 = 104000 + np.array([205, 208, 209, 210, 211, 212])
	anode_strt2 = np.array([7, 3, 7, 5, 7, 7])
	anode_duration2 = np.array([2, 6, 4, 4, 6, 8])
	
	# GVA=20psi, GVC=60psi, GVN=100psi
	shots_anode3 = 104000 + np.array([237, 238, 239, 240, 241, 242])
	anode_strt3 = [7, 5, 3, 7, 7, 7]
	anode_duration3 = [2, 4, 6, 4, 6, 8]
	
	# GVA=40psi, GVC=60psi, GVN=100psi
	shots_anode4 = 104000 + np.array([243, 245, 246, 251, 249, 250])
	anode_strt4 = [7, 5, 3, 7, 7, 7]
	anode_duration4 = [2, 4, 6, 4, 6, 8]
	
	# GVA=25psi, GVC=60psi, GVN=100psi
	shots_anode5 = 104000 + np.array([253, 254, 255, 252, 256, 257])
	anode_strt5 = [7, 5, 3, 7, 7, 7]
	anode_duration5 = [2, 4, 6, 4, 6, 8]
	
	# GVA=30psi, GVC=60psi, GVN=100psi
	shots_cathode = 104000 + np.array([177, 193, 195, 196, 197, 200, 201, 202, 203, 204])
	cathode_strt = [6, 6, 4, 2, 0, 8, 8, 6, 6, 6]
	cathode_duration = [4, 4, 6, 8, 10, 2, 2, 2, 7, 9]
	shots_cathode_noneut = 104000 + np.array([199])
	cathode_strt_noneut = np.array([0])
	cathode_duration_noneut = np.array([10])
	
	# GVA=30psi, GVC=80psi, GVN=100psi
	shots_cathode2 = 104000 + np.array([213, 214, 215, 216, 222])
	cathode_strt2 = [6, 4, 2, 6, 6]
	cathode_duration2 = [4, 6, 8, 2, 6]
	
	# GVA=30psi, GVC=40psi, GVN=100psi
	shots_cathode3 = 104000 + np.array([224, 226, 227, 228, 230])
	cathode_strt3 = [6, 4, 2, 6, 6]
	cathode_duration3 = [4, 6, 8, 2, 6]
	
	# GVA=30psi, GVC=20psi, GVN=100psi
	shots_cathode4 = 104000 + np.array([258, 259, 260, 261, 262])
	cathode_strt4 = [6, 4, 2, 6, 6]
	cathode_duration4 = [4, 6, 8, 2, 6]
	
	# GVA=30psi, GVC=50psi, GVN=100psi
	shots_cathode5 = 104000 + np.array([263, 264, 265, 266, 268])
	cathode_strt5 = [6, 4, 2, 6, 6]
	cathode_duration5 = [4, 6, 8, 2, 6]
	
	fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True)
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	plot_some_data(anode_duration3, shots_anode3, clrs[0], '20/60/100', ax1)
	plot_some_data(anode_duration5, shots_anode5, clrs[1], '25/60/100', ax1)
	plot_some_data(anode_duration, shots_anode, clrs[2], '30/60/100', ax1)
	plot_some_data(anode_duration4, shots_anode4, clrs[3], '40/60/100', ax1)
	plot_some_data(anode_duration2, shots_anode2, clrs[4], '60/60/100', ax1)
	# plot_some_data(anode_duration_noneut, shots_anode_noneut, clrs[0], 'anode no neutralizer')
	plot_some_data(cathode_duration4, shots_cathode4, clrs[0], '30/20/100', ax2)
	plot_some_data(cathode_duration3, shots_cathode3, clrs[1], '30/40/100', ax2)
	plot_some_data(cathode_duration5, shots_cathode5, clrs[2], '30/50/100', ax2)
	plot_some_data(cathode_duration, shots_cathode, clrs[3], '30/60/100', ax2)
	# plot_some_data(cathode_duration_noneut, shots_cathode_noneut, clrs[1], 'cathode no neutralizer')
	plot_some_data(cathode_duration2, shots_cathode2, clrs[4], '30/80/100', ax2)
	
	# ax1.axvline(5, ls='--', c='k', alpha=0.5)
	# ax1.annotate('nominal', (5, .5))
	ax1.plot([0], [0])
	ax1.spines.left.set_position('zero')
	ax1.spines.right.set_color('none')
	ax1.spines.bottom.set_position('zero')
	ax1.spines.top.set_color('none')
	ax1.xaxis.set_ticks_position('bottom')
	ax1.yaxis.set_ticks_position('left')
	ax1.set_xlabel('puff duration (ms)')
	ax1.set_ylabel('measured/predicted dT')
	ax1.set_title('Anode')
	ax1.legend(title='GVA/GVC/GVN', fontsize=10)
	ax2.plot([0], [0])
	ax2.spines.left.set_position('zero')
	ax2.spines.right.set_color('none')
	ax2.spines.bottom.set_position('zero')
	ax2.spines.top.set_color('none')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.yaxis.set_ticks_position('left')
	ax2.set_xlabel('puff duration (ms)')
	ax2.set_title('Cathode')
	ax2.legend(title='GVA/GVC/GVN', fontsize=10)


def new_valve_data():
	shots7Dec = 104500 + np.array(
		[10, 12, 13, 14, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 37, 40, 42, 43, 44, 45, 46])
	hv = np.array(
		[15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 10, 10, 10, 10, 10, 10])
	i15kv = np.where(hv == 15)[0]
	i17kv = np.where(hv == 17)[0]
	i10kv = np.where(hv == 10)[0]
	dt_pred, dt_meas = [], []  # shot arrays for joules
	perv, dperv = [], []
	for ish, shot in enumerate(shots7Dec):
		dtm, dtp = cal_dtemp(shot)
		dt_meas.append(dtm)
		dt_pred.append(dtp)
		p, dp, _, _ = calculate_perv_pwr(shot)
		perv.append(p)
		dperv.append(dp)
	
	fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(14, 5))
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	slope = np.nanmean([dt_meas[i] / dt_pred[i] for i in range(len(dt_meas))])
	
	dt_pred, dt_meas, perv, dperv = np.array(dt_pred), np.array(dt_meas), np.array(perv), np.array(dperv)
	ax1.plot([0, max(dt_pred)], [0, max(dt_pred) * slope], '--', c=clrs[3])
	ax1.plot(dt_pred[i15kv], dt_meas[i15kv], 'o', label='HV=15kV', c=clrs[0])
	ax1.plot(dt_pred[i17kv], dt_meas[i17kv], 'o', label='HV=17kV', c=clrs[1])
	ax1.plot(dt_pred[i10kv], dt_meas[i10kv], 'o', label='HV=10kV', c=clrs[2])
	
	ax1.annotate(f'{slope * 100.:.2f}%', (max(dt_pred), max(dt_pred) * slope), c=clrs[0])
	
	ax1.set_title('20psi, cathode only 6-10 ms')
	ax1.set_xlabel('predicted $\Delta T$ ($\degree$C)')
	ax1.set_ylabel('measured $\Delta T$ ($\degree$C)')
	ax1.set_title('Calorimeter')
	ax1.legend(fontsize=10)
	
	ax2.errorbar(perv[i15kv] * 1.e6, dt_meas[i15kv] / dt_pred[i15kv], xerr=dperv[i15kv] * 1.e6, marker='o', ls='',
	             label='HV=15kV', c=clrs[0])
	ax2.errorbar(perv[i17kv] * 1.e6, dt_meas[i17kv] / dt_pred[i17kv], xerr=dperv[i17kv] * 1.e6, marker='o', ls='',
	             label='HV=17kV', c=clrs[1])
	ax2.errorbar(perv[i10kv] * 1.e6, dt_meas[i10kv] / dt_pred[i10kv], xerr=dperv[i10kv] * 1.e6, marker='o', ls='',
	             label='HV=10kV', c=clrs[2])
	ax2.set_xlabel('perv (e-6)')
	ax2.set_ylabel('measured fraction of predicted $\Delta T$')


def new_valve_neutralizer():
	shots7Dec = 104500 + np.array([66, 68, 69, 70, 72, 73, 74, 75, 76, 77, 78, 80, 81, 82])
	neut_puffstart = np.array([15, 15, 14, 14, 14, 14, 13, 13, 12, 12, 11, 11, 10, 10])
	neut_psi = np.array([15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15])
	dt_pred, dt_meas = [], []  # shot arrays for joules
	perv, dperv = [], []
	for ish, shot in enumerate(shots7Dec):
		dtm, dtp = cal_dtemp(shot)
		dt_meas.append(dtm)
		dt_pred.append(dtp)
		p, dp, _, _ = calculate_perv_pwr(shot)
		perv.append(p)
		dperv.append(dp)
	
	fig, (ax1, ax2) = plt.subplots(ncols=2)  # , figsize=(12, 6))
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	slope = np.nanmean([dt_meas[i] / dt_pred[i] for i in range(len(dt_meas))])
	
	psi15 = np.where(neut_psi == 15)
	dt_pred, dt_meas, perv, dperv = np.array(dt_pred), np.array(dt_meas), np.array(perv), np.array(dperv)
	ax1.plot(neut_puffstart[psi15], dt_meas[psi15] / dt_pred[psi15], 'o', label='psi=15', c=clrs[0])
	ax1.annotate(f'{slope * 100.:.2f}%', (max(dt_pred), max(dt_pred) * slope), c=clrs[0])
	
	ax1.set_title('15psi, variable neutralizer puff start')
	ax1.set_xlabel('neutralizer puff start (ms)')
	ax1.set_ylabel('measured fraction of predicted $\Delta T$')
	ax1.set_title('Calorimeter')
	ax1.legend(fontsize=10)
	ax1.set_xlim(left=0)
	ax1.set_ylim(bottom=0)
	
	ax2.errorbar(perv[psi15] * 1.e6, dt_meas[psi15] / dt_pred[psi15], xerr=dperv[psi15] * 1.e6, marker='o', ls='',
	             c=clrs[0])
	ax2.set_xlabel('perv (e-6)')
	ax2.set_ylabel('measured fraction of predicted $\Delta T$')
	ax2.set_xlim(left=0)
	ax2.set_ylim(bottom=0)


def new_valve_cathode():
	# ADD PLOT OF Ibeam SLOPE VS SKEW to show "shape" of current traces vs parameters
	cath_scan_shots = 104000 + np.array(
		[584, 585, 586, 587, 589, 598, 599, 613, 617, 624, 628, 629, 630, 631, 632, 643, 642, 641, 640, 638, 637, 636,
		 649, 648, 646, 650, 651, 660, 661, 662, 663, 664, 665, 674, 673, 672, 671, 670, 669, 668, 713, 714, 712, 710,
		 709, 703, 701, 699, 707, 718, 717, 716, 715])
	puffstart = np.array(
		[4, 4, 4, 4, 4, 6, 6, 6, 6, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0,
		 0, 0, 0, 6, 6, 4, 4, 4, 2, 2, 2, 2, 0, 0, 0, 0])
	puffduration = np.array(
		[6, 4, 2, 8, 10, 8, 6, 4, 2, 8, 6, 4, 2, 10, 12, 2, 4, 6, 8, 10, 12, 14, 2, 4, 6, 8, 10, 2, 4, 6, 8, 10, 12, 2,
		 4, 6, 8, 10, 12, 14, 4, 8, 2, 6, 10, 2, 6, 8, 12, 2, 6, 10, 14])
	psi = np.array(
		[20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 30, 30, 30, 30, 30, 30,
		 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10])
	dt_pred, dt_meas = [], []  # shot arrays for joules
	perv, dperv, pwr, dpwr, dt0, dt1, dt2, dt3, dt4 = [], [], [], [], [], [], [], [], []
	ipslope, ipstdev = [], []
	for ish, shot in enumerate(cath_scan_shots):
		dtm, dtp, dt_arr, ip_shape = cal_dtemp(shot, more=True)
		dt0.append(dt_arr[0])
		dt1.append(dt_arr[1])
		dt2.append(dt_arr[2])
		dt3.append(dt_arr[3])
		dt4.append(dt_arr[4])
		ipslope.append(ip_shape[0])
		ipstdev.append(ip_shape[1])
		dt_meas.append(dtm)
		dt_pred.append(dtp)
		p, dp, pow, dpow = calculate_perv_pwr(shot)
		perv.append(p)
		dperv.append(dp)
		pwr.append(pow)
		dpwr.append(dpow)
	
	fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(12, 5))
	fig2, (axa, axb, axc, axd, axe) = plt.subplots(ncols=5, figsize=(12, 5), sharey=True)
	fig3, axz = plt.subplots()
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	slope = np.nanmean([dt_meas[i] / dt_pred[i] for i in range(len(dt_meas))])
	
	dt_pred, dt_meas, perv, dperv = np.array(dt_pred), np.array(dt_meas), np.array(perv), np.array(dperv)
	dt0, dt1, dt2, dt3, dt4 = np.array(dt0), np.array(dt1), np.array(dt2), np.array(dt3), np.array(dt4)
	pwr, dpwr, ipslope, ipstdev = np.array(pwr), np.array(dpwr), np.array(ipslope), np.array(ipstdev)
	for ps in np.unique(puffstart):
		if ps == 0:
			mrk = 'x'
		if ps == 2:
			mrk = 'd'
		if ps == 4:
			mrk = 'o'
		if ps == 6:
			mrk = 's'
		
		psi20 = np.where((psi == 20) & (puffstart == ps))
		ax1.plot(puffduration[psi20], dt_meas[psi20] / dt_pred[psi20], mrk, c=clrs[0])
		ax2.errorbar(pwr[psi20] * 1.e-3, perv[psi20] * 1.e6, xerr=dpwr[psi20] * 1.e-3, yerr=dperv[psi20] * 1.e6,
		             marker=mrk, ls='', c=clrs[0])
		ax3.plot(ipstdev[psi20], ipslope[psi20], mrk, c=clrs[0])
		for i, (ax, thermocoup) in enumerate(zip([axa, axb, axc, axd, axe], [dt0, dt1, dt2, dt3, dt4])):
			for pdval, dtval in zip(puffduration[psi20], thermocoup[psi20]):
				ax.plot(pdval, dtval, mrk, c=clrs[0])
		
		psi30 = np.where((psi == 30) & (puffstart == ps))
		ax1.plot(puffduration[psi30], dt_meas[psi30] / dt_pred[psi30], mrk, c=clrs[2])
		ax2.errorbar(pwr[psi30] * 1.e-3, perv[psi30] * 1.e6, xerr=dpwr[psi30] * 1.e-3, yerr=dperv[psi30] * 1.e6,
		             marker=mrk, ls='', c=clrs[2])
		ax3.plot(ipstdev[psi30], ipslope[psi30], mrk, c=clrs[2])
		for i, (ax, thermocoup) in enumerate(zip([axa, axb, axc, axd, axe], [dt0, dt1, dt2, dt3, dt4])):
			for pdval, dtval in zip(puffduration[psi30], thermocoup[psi30]):
				ax.plot(pdval, dtval, mrk, c=clrs[2])
		
		psi10 = np.where((psi == 10) & (puffstart == ps))
		ax1.plot(puffduration[psi10], dt_meas[psi10] / dt_pred[psi10], mrk, c=clrs[1])
		ax2.errorbar(pwr[psi10] * 1.e-3, perv[psi10] * 1.e6, xerr=dpwr[psi10] * 1.e-3, yerr=dperv[psi10] * 1.e6,
		             marker=mrk, ls='', c=clrs[1])
		ax3.plot(ipstdev[psi10], ipslope[psi10], mrk, c=clrs[1])
		for i, (ax, thermocoup) in enumerate(zip([axa, axb, axc, axd, axe], [dt0, dt1, dt2, dt3, dt4])):
			for pdval, dtval in zip(puffduration[psi10], thermocoup[psi10]):
				ax.plot(pdval, dtval, mrk, c=clrs[1])
		
		axz.plot(ipslope[psi10], dt0[psi10] / dt1[psi10], mrk, c=clrs[1])
		axz.plot(ipslope[psi20], dt0[psi20] / dt1[psi20], mrk, c=clrs[0])
		axz.plot(ipslope[psi30], dt0[psi30] / dt1[psi30], mrk, c=clrs[2])
		if ps == 4:
			print(f'avg center/right ratio for PSI=20, PS=4 is {np.mean(dt0[psi20] / dt1[psi20]):.2f}')
	
	ax1.set_title('cathode puff timing/duration scan')
	ax1.set_xlabel('neutralizer puff duration (ms)')
	ax1.set_ylabel('measured fraction of predicted $\Delta T$')
	ax1.set_title('Calorimeter')
	ax1.set_xlim(left=0)
	ax1.set_ylim((0, .6))
	m1, = ax1.plot(np.nan, np.nan, 'kx')
	m2, = ax1.plot(np.nan, np.nan, 'kd')
	m3, = ax1.plot(np.nan, np.nan, 'ko')
	m4, = ax1.plot(np.nan, np.nan, 'ks')
	leg1 = ax1.legend([m1, m2, m3, m4], ['0', '2', '4', '6'], fontsize=10, title='puff start time', loc='lower right')
	p1, = ax1.plot(np.nan, np.nan, '-', c=clrs[0], label='20')
	p2, = ax1.plot(np.nan, np.nan, '-', c=clrs[1], label='10')
	p3, = ax1.plot(np.nan, np.nan, '-', c=clrs[2], label='30')
	ax1.legend(handles=[p1, p2, p3], title='psi', loc='lower center')
	ax1.add_artist(leg1)
	
	ax2.axhline(15, ls='--', c='k', label='target')
	ax2.legend()
	ax2.set_xlabel('beam power (kW)')
	ax2.set_ylabel('perv (e-6)')
	ax2.set_xlim(left=0)
	ax2.set_ylim(bottom=0)
	
	ax3.set_xlabel('stdev')
	ax3.set_ylabel('slope of linear fit')
	
	axa.set_ylabel('$\Delta T$')
	axc.set_xlabel('puff duration (ms)')
	axa.set_title('Center')
	axb.set_title('Right')
	axc.set_title('Up')
	axd.set_title('Left')
	axe.set_title('Down')
	
	axz.set_xlabel('Slope of Ibeam')
	axz.set_ylabel('Center/Right temp increase ratio')
	
	plt.tight_layout()


def calorimeter_21Jan():
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	
	def do_stuff(shots, mark, inlegend=True):
		dt_pred, dt_meas = [], []
		perv, dperv, pwr, dpwr, dt0, dt1, dt2, dt3, dt4 = [], [], [], [], [], [], [], [], []
		for ish, shot in enumerate(shots):
			dtm, dtp, dt_arr, _ = cal_dtemp(shot, more=True)
			avp = avg_perv(shot)
			dt0.append(dt_arr[0])
			dt1.append(dt_arr[1])
			dt2.append(dt_arr[2])
			dt3.append(dt_arr[3])
			dt4.append(dt_arr[4])
			dt_meas.append(dtm)
			dt_pred.append(dtp)
			perv.append(avp)
		dt_pred, dt_meas, perv, dperv = np.array(dt_pred), np.array(dt_meas), np.array(perv), np.array(dperv)
		print(
			f'max: {np.max(dt_meas / dt_pred)} on shot {shots[np.where(dt_meas / dt_pred == np.max(dt_meas / dt_pred))[0][0]]}')
		dt0, dt1, dt2, dt3, dt4 = np.array(dt0), np.array(dt1), np.array(dt2), np.array(dt3), np.array(dt4)
		# pwr, dpwr, ipslope, ipstdev = np.array(pwr), np.array(dpwr), np.array(ipslope), np.array(ipstdev)
		perv *= 1.e6
		
		xx = np.linspace(min(perv), max(perv))
		
		ax1.plot(perv, dt_meas / dt_pred, mark, c=clrs[0])
		if inlegend:
			lbls = ['center', 'right', 'up', 'left', 'down']
		else:
			lbls = ['_nolegend_'] * 5
		ax2.plot(perv, dt0, mark, label=lbls[0], c=clrs[0])
		ax2.plot(perv, dt1, mark, label=lbls[1], c=clrs[1])
		ax2.plot(perv, dt2, mark, label=lbls[2], c=clrs[2])
		ax2.plot(perv, dt3, mark, label=lbls[3], c=clrs[3])
		ax2.plot(perv, dt4, mark, label=lbls[4], c=clrs[4])
		ax3.plot(perv, dt0 / dt1, mark, c=clrs[0])
		print(f'max power fraction on calorimeter {max(dt_meas / dt_pred) * 100}%')
	
	# NEW ORIFICES 01Dec21
	fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(12, 5))
	shots21jan22 = 104800 + np.array([66, 68, 69, 70, 71, 72, 74, 75, 76, 77, 78, 80, 85, 86, 87, 89, 92, 93])
	do_stuff(shots21jan22, 'o')
	shots9dec21 = 104800 + np.array([97, 98, 99])  # replicating shots 104584-89 to compare ratio between center/right
	# do_stuff(shots9dec21, 's', inlegend=False)
	shots25jan22 = 104900 + np.array([26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 43, 44, 45, 46])
	do_stuff(shots25jan22, 'd', inlegend=False)
	shots25jan22_noneut = 104900 + np.array([38, 41, 48, 49, 50, 51])
	# do_stuff(shots25jan22_noneut, 's')
	ax1.set_ylabel('Fractional power onto cal')
	ax2.set_xlabel('perveance (e-6)')
	ax2.set_ylabel('dT')
	ax2.legend()
	ax3.set_ylabel('dT ratio center/right')
	for ax in [ax1, ax2, ax3]:
		ax.set_ylim(bottom=0)
	
	plt.tight_layout()


plt.show()


def get_thermocouple_numbers(shots, tp, ret_minmax=False):
	c0p, c1p, c2p, c3p, c4p = np.zeros_like(tp), np.zeros_like(tp), np.zeros_like(tp), np.zeros_like(
		tp), np.zeros_like(tp)
	num = 0.
	minmax = np.zeros((5, 2))
	minmax[:] = np.nan
	for shot in shots:
		if 12.5 < avg_perv(shot) * 1.e6 < 18:
			num += 1
			t, c0, c1, c2, c3, c4 = cal_temp_sigs(shot)
			c0p += np.interp(tp, t, c0)
			c1p += np.interp(tp, t, c1)
			c2p += np.interp(tp, t, c2)
			c3p += np.interp(tp, t, c3)
			c4p += np.interp(tp, t, c4)
			ipeak = np.where(c0 == max(c0))[0][0]
			for i in np.arange(5):
				minmax[i, 0] = np.nanmin([minmax[i, 0], [c0, c1, c2, c3, c4][i][ipeak]])
				minmax[i, 1] = np.nanmax([minmax[i, 1], [c0, c1, c2, c3, c4][i][ipeak]])
		else:
			print(f'bad perveance for shot {shot}')
	c0p, c1p, c2p, c3p, c4p = c0p / num, c1p / num, c2p / num, c3p / num, c4p / num
	ipeak = np.where(c0p == max(c0p))[0][0]
	# dt0, dt1, dt2, dt3, dt4 = max(c0p), max(c1p), max(c2p), max(c3p), max(c4p)
	dt0, dt1, dt2, dt3, dt4 = c0p[ipeak], c1p[ipeak], c2p[ipeak], c3p[ipeak], c4p[ipeak]
	if ret_minmax:
		return c0p, c1p, c2p, c3p, c4p, [dt0, dt1, dt2, dt3, dt4], minmax
	else:
		return c0p, c1p, c2p, c3p, c4p, [dt0, dt1, dt2, dt3, dt4]


def compare_neutralizer_onoff():
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	tp = np.linspace(-50, 100, endpoint=True)
	# new orifice, realignment 1, shift to right (aim left)
	shots25jan22 = 104900 + np.array([26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 43, 44, 45, 46])
	# new orifice, new alignment, no neutralizer
	shots25jan22_noneut = 104900 + np.array([38, 41, 48, 49, 50, 51])
	pre0, pre1, pre2, pre3, pre4, dt_pre = get_thermocouple_numbers(shots25jan22, tp)
	post0, post1, post2, post3, post4, dt_post = get_thermocouple_numbers(shots25jan22_noneut, tp)
	fig = plt.figure(tight_layout=True, figsize=(10, 8))
	gs = gridspec.GridSpec(2, 2)
	ax1 = fig.add_subplot(gs[0, 0])
	ax2 = fig.add_subplot(gs[0, 1])
	ax3 = fig.add_subplot(gs[1, :])
	ax1.plot(tp, pre0, clrs[0], label='center')  # correct labeling
	ax1.plot(tp, pre1, clrs[1], label='down')
	ax1.plot(tp, pre2, clrs[2], label='right')
	ax1.plot(tp, pre3, clrs[3], label='up')
	ax1.plot(tp, pre4, clrs[4], label='left')
	ax2.plot(tp, post0, clrs[0])
	ax2.plot(tp, post1, clrs[1])
	ax2.plot(tp, post2, clrs[2])
	ax2.plot(tp, post3, clrs[3])
	ax2.plot(tp, post4, clrs[4])
	ax1.legend()
	
	for i in range(5):
		ax3.plot([-1, 1], [dt_pre[i], dt_post[i]], 's-', c=clrs[i])
	
	ax1.set_ylabel('thermocouple\ntemp (deg C)')
	ax1.set_xlabel('time rel to beam')
	ax2.set_ylim(ax1.get_ylim())
	ax2.set_xlabel('time rel to beam')
	ax3.set_ylabel('thermocouple\n$\Delta T$ (deg C)')
	ax3.set_xlim((-2, 2))
	ax3.set_xticks([-1, 1])
	
	ax1.set_title('Neutralizer ON')
	ax2.set_title('Neutralizer OFF')
	ax3.set_xticklabels(['neut ON', 'neut OFF'])
	
	plt.show()


def beam_realignment(all4=False):
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	
	# new orifice, original alignment
	shots21jan22 = 104800 + np.array([66, 68, 69, 70, 71, 72, 74, 75, 76, 77, 78, 80, 85, 86, 87, 89, 92, 93])
	# new orifice, realignment 1, shift to right (aim left)
	shots25jan22 = 104900 + np.array([26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 43, 44, 45, 46])
	# new orifice, realignment 2, shift part way back left (aim right)
	shots1feb22 = 105000 + np.array([9, 11, 12, 13, 14, 15])
	# realignment 3: centered calorimeter based on measurements made on spare
	shots3feb22 = 105000 + np.array([22, 23, 24, 25, 27])
	tp = np.linspace(-50, 100, endpoint=True)
	
	a0, a1, a2, a3, a4, adt = get_thermocouple_numbers(shots21jan22, tp)
	if all4:
		b0, b1, b2, b3, b4, bdt = get_thermocouple_numbers(shots25jan22, tp)
		c0, c1, c2, c3, c4, cdt = get_thermocouple_numbers(shots1feb22, tp)
	d0, d1, d2, d3, d4, ddt = get_thermocouple_numbers(shots3feb22, tp)
	
	if all4:
		fig = plt.figure(tight_layout=True, figsize=(10, 8))
		gs = gridspec.GridSpec(2, 4)
		ax1 = fig.add_subplot(gs[0, 0])
		ax2 = fig.add_subplot(gs[0, 1])
		ax3 = fig.add_subplot(gs[0, 2])
		ax4 = fig.add_subplot(gs[0, 3])
		bax = fig.add_subplot(gs[1, :])
	else:
		fig = plt.figure(figsize=(6, 5))
		gs = gridspec.GridSpec(2, 2)
		ax1 = fig.add_subplot(gs[0, 0])
		ax4 = fig.add_subplot(gs[0, 1])
		bax = fig.add_subplot(gs[1, :])
	lbls = ['center', 'down', 'right', 'up', 'left']  # correct labeling
	for i, a in enumerate([a0, a1, a2, a3, a4]):
		ax1.plot(tp, a, clrs[i], label=lbls[i])
	if all4:
		for i, b in enumerate([b0, b1, b2, b3, b4]):
			ax2.plot(tp, b, clrs[i])
		for i, c in enumerate([c0, c1, c2, c3, c4]):
			ax3.plot(tp, c, clrs[i])
	for i, d in enumerate([d0, d1, d2, d3, d4]):
		ax4.plot(tp, d, clrs[i])
	
	for i in range(5):
		if all4:
			xt = [0, 1, 2, 3]
			bax.plot(xt, [adt[i], bdt[i], cdt[i], ddt[i]], 's-', c=clrs[i])
		else:
			xt = [0, 1]
			bax.plot([0, 1], [adt[i], ddt[i]], 's-', c=clrs[i], label=lbls[i])
	
	fs, fs2 = 12, 10
	ax1.set_ylabel('thermocouple\ntemp (deg C)', fontsize=fs)
	ax1.set_xlabel('time rel to beam', fontsize=fs)
	if all4:
		for ax in [ax2, ax3, ax4]:
			ax.set_ylim(ax1.get_ylim())
	else:
		ax4.set_ylim(ax1.get_ylim())
		for ax in [ax1, ax4, bax]:
			ax.tick_params(labelsize=fs)
	# ax2.set_xlabel('time rel to beam')
	bax.set_ylabel('thermocouple\n$\Delta T$ (deg C)', fontsize=fs)
	bax.set_xlim((min(xt) - .5, max(xt) + .5))
	bax.set_xticks(xt)
	
	# COMPUTE rough adjustment guess based on thermocouple ratios
	if all4:
		fig2, axx = plt.subplots()
		xx = [0, .1, .3]  # c, a, b
		yy = [abs(cdt[4] - cdt[2]), abs(adt[4] - adt[2]),
		      abs(bdt[4] - bdt[2])]  # left/right asymmetry across 2 adjustments
		yy2 = abs(ddt[3] - ddt[1])  # up/down asymmetry
		fit = np.polyfit(xx, yy, 1)
		axx.plot(xx, yy, 'o-')
		axx.plot((yy2 - fit[1]) / fit[0], yy2, 'rs')
		print(f'Suggested move {yy2}"')
	# lr_pre, lr_post = dt_pre[4] - dt_pre[2], dt_post[4] - dt_post[2]
	# crossing_fraction = abs(lr_pre) / (abs(lr_pre) + abs(lr_post))  # fraction between pre/post where L/R cross
	# print(f'L/R cross at {crossing_fraction * 100:.2f}% between Before and After')
	
	ax1.set_title('Original', fontsize=fs)
	if all4:
		ax2.set_title('Source adjusted 0.4" right')
		ax3.set_title('Source 0.3" left')
		ax4.set_title('Calorimeter lowered ~0.4"')
	else:
		ax4.set_title('Centered', fontsize=fs)
	bax.set_xlabel('beam realignment', fontsize=fs)
	if all4:
		bax.set_xticklabels(['original', 'adj1', 'adj2', 'adj3'])
	else:
		bax.set_xticklabels(['original', 'centered'], fontsize=fs2)
	bax.legend(title='thermocouple', loc='right', fontsize=fs2, frameon=False)
	plt.tight_layout()
	plt.show()


def calorimeter_positional_variation():
	import matplotlib.gridspec as gridspec
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	
	tp = np.linspace(-50, 100, endpoint=True)
	
	def get_numbers(shots):
		therms = np.zeros((len(tp), 5))
		num = 0.
		for shot in shots:
			num += 1
			t, c0, c1, c2, c3, c4 = cal_temp_sigs(shot, calonly=True)
			therms[:, 0] += np.interp(tp, t, c0)
			therms[:, 1] += np.interp(tp, t, c1)
			therms[:, 2] += np.interp(tp, t, c2)
			therms[:, 3] += np.interp(tp, t, c3)
			therms[:, 4] += np.interp(tp, t, c4)
		therms /= num
		ipeak = np.where(therms[:, 0] == max(therms[:, 0]))[0][0]
		dt = np.zeros(5)
		for i in range(5):
			dt[i] = therms[ipeak, i]
		return therms, dt
	
	shots1 = 104900 + np.array([75, 76, 77, 78, 80])  # nominal- calorimeter centered on beamline
	shots2 = 104900 + np.array([81, 82, 83, 84, 85])  # position +0.25"
	shots3 = 104900 + np.array([86, 88, 89, 90, 92])  # position at +0.5"
	shots4 = 104900 + np.array([96, 97, 98, 99, 101])
	shots5 = 105000 + np.array([2, 3, 4, 5, 6])
	p1, dt1 = get_numbers(shots1)
	p2, dt2 = get_numbers(shots2)
	p3, dt3 = get_numbers(shots3)
	p4, dt4 = get_numbers(shots4)
	p5, dt5 = get_numbers(shots5)
	fig = plt.figure(tight_layout=True)  # , figsize=(10, 8))
	gs = gridspec.GridSpec(2, 5)
	ax1 = fig.add_subplot(gs[0, 0])
	ax2 = fig.add_subplot(gs[0, 1])
	ax3 = fig.add_subplot(gs[0, 2])
	ax4 = fig.add_subplot(gs[0, 3])
	ax5 = fig.add_subplot(gs[0, 4])
	axb = fig.add_subplot(gs[1, :])
	lbls = ['center', 'down', 'right', 'up', 'left']
	for ip in range(5):
		ax1.plot(tp, p1[:, ip], clrs[ip])
		ax2.plot(tp, p2[:, ip], clrs[ip])
		ax3.plot(tp, p3[:, ip], clrs[ip])
		ax4.plot(tp, p4[:, ip], clrs[ip])
		ax5.plot(tp, p5[:, ip], clrs[ip], label=f'{ip}: {lbls[ip]}')
	ax5.legend()
	
	pos = [0, .25, .5, .75, 1.0]
	axb.set_xticks(pos)
	axb.set_xticklabels([str(p) for p in pos])
	for i in range(5):
		axb.plot(pos, [dt1[i], dt2[i], dt3[i], dt4[i], dt5[i]], 's-', c=clrs[i])
	
	ax1.set_ylabel('thermocouple\ntemp (deg C)')
	ax3.set_xlabel('time rel to beam')
	for ax in [ax1, ax2, ax4, ax5]:
		ax.set_ylim(ax3.get_ylim())
	axb.set_ylabel('thermocouple\n$\Delta T$ (deg C)')
	
	for ia, ax in enumerate([ax1, ax2, ax3, ax4, ax5]):
		ax.set_title(f'{pos[ia]}')
	axb.set_xlabel('calorimeter position')
	
	plt.show()


def perveance_scan_7feb22(plot_vs='fwhm'):
	"""
	beam aligned on calorimeter, plot fwhm vs perveance
	plot_vs: fwhm, power_frac_to_cal, or efold
	"""
	if plot_vs.lower() == 'fwhm' or plot_vs.lower() == 'efold':
		units = ' (mm)'
	else:
		units = ''
	
	sh1 = 105100 + np.array(
		[24, 27, 31, 32, 35, 45, 47, 48, 49, 50, 51, 52, 53, 54, 55, 57, 58, 59, 61, 62, 63, 64, 65, 67, 68, 69, 71, 73,
		 74, 75, 76])
	# sh1 = 105100 + np.array([45, 47, 50, 55, 58, 61, 62, 63, 67, 68, 73, 74, 75, 76])
	fig, (ax1) = plt.subplots(ncols=1)  # 2, figsize=(10, 5))
	fwhm_arr, perv_arr = np.array([]), np.array([])
	insuff = []
	for ish, shot in enumerate(sh1):
		perv_arr = np.append(perv_arr, avg_perv(shot))
		fwhm, suff_neut = cal_gauss_fit(shot, ret=plot_vs)
		if not suff_neut:
			insuff.append(ish)
		fwhm_arr = np.append(fwhm_arr, fwhm)
	
	fs = 14
	ax1.plot(perv_arr * 1.e6, fwhm_arr, 'o')
	# ax1.plot(perv_arr[insuff] * 1.e6, fwhm_arr[insuff], 'kx')
	ax1.set_xlabel('perveance ($\\times 10^{-6}$)', fontsize=fs)
	ax1.set_ylabel(f'{plot_vs}{units}', fontsize=fs)
	ax1.tick_params(labelsize=fs)
	ax1.set_xlim((0, 30))
	plt.tight_layout()


def special_calibration(redo=False):
	lvm_files = ['Y:/thermocouple/Calorimeter/04082022171755.lvm',
	             'Y:/thermocouple/Calorimeter/04062022125726.lvm',
	             'Y:/thermocouple/Calorimeter/04062022110632.lvm',
	             'Y:/thermocouple/Calorimeter/04012022085324.lvm',
	             'Y:/thermocouple/Calorimeter/03292022132144.lvm',
	             'Y:/thermocouple/Calorimeter/03292022141515.lvm',
	             'Y:/thermocouple/Calorimeter/04012022094257.lvm']
	calibs = [-27, -20, -43, 97, 26, 26, 97]
	if redo:
		# NEED DIFFERENT CALIBRATIONS for each lvm file <sigh>
		# direc = 'Y:/thermocouple/Calorimeter/'
		# lvm_files = os.listdir(direc)
		# lvm_files = [f'{direc}{lvmf}' for lvmf in lvm_files if '2022' in lvmf]  # lvms with date format
		shots = np.arange(508397, 508398)
		for ifile, lvmf in enumerate(lvm_files):
			print(f'analyzing {lvmf}')
			shotsinfile = []
			shot_offsets = []
			for shot in shots:
				timestamp = get_shot_timestamp(shot)
				# print(f'shot {shot} occurred at {timestamp}')
				shot_ts = datetime.datetime.strptime(timestamp, '%m/%d/%Y %I:%M:%S %p')
				lvm = lvm_read.read(lvmf, dump_file=False, read_from_pickle=False)
				# lvm_strt = datetime.datetime.strptime(lvm['Date'] + ':' + lvm['Time'].split('.')[0], '%Y/%m/%d:%H:%M:%S')  # CLOCK IS OFF???
				lvm_strt = datetime.datetime.strptime(lvmf.split('/')[-1].split('.')[0], '%m%d%Y%H%M%S')
				lvm0 = lvm[0]
				lvmtime, c0 = lvm0['data'][:, 0], lvm0['data'][:, 15]
				lvm_end = lvm_strt + datetime.timedelta(seconds=lvmtime[-1])
				if lvm_strt < shot_ts < lvm_end:
					shot_offsets.append((shot_ts - lvm_strt).seconds)
					shotsinfile.append(shot)
			if len(shot_offsets) > 0:
				plt.plot(lvmtime, c0)
				for shoff in shot_offsets:
					ii = np.where((lvmtime + calibs[ifile] > shoff - 1) & (shoff + 1 > lvmtime + calibs[ifile]))
					plt.plot(lvmtime[ii], c0[ii], 'r')
				a = 1
			else:
				print(f'no shots found in {lvmf}')
	return lvm_files, calibs


def special_extraction(lvm_files, shot, calibs):
	# todo: NEED TO VERIFY THIS IS PROPERLY IDENTIFYING SIGNAL INTERVALS
	timestamp = get_shot_timestamp(shot)
	shot_ts = datetime.datetime.strptime(timestamp, '%m/%d/%Y %I:%M:%S %p')
	for (cal, lvf) in zip(calibs, lvm_files):
		lvm = lvm_read.read(lvf, dump_file=False, read_from_pickle=False)
		# lvm_strt = datetime.datetime.strptime(lvm['Date'] + ':' + lvm['Time'].split('.')[0], '%Y/%m/%d:%H:%M:%S')  # CLOCK IS OFF???
		lvm_strt = datetime.datetime.strptime(lvf.split('/')[-1].split('.')[0], '%m%d%Y%H%M%S')
		lvm0 = lvm[0]
		lvmtime = lvm0['data'][:, 0]
		lvm_end = lvm_strt + datetime.timedelta(seconds=lvmtime[-1])
		if lvm_strt < shot_ts < lvm_end:
			# print(f'found shot {shot} in lvm file {lvf}')
			t0 = (shot_ts - lvm_strt).seconds
			ii = np.where((lvmtime + cal > t0 - 20) & (t0 + 90 > lvmtime + cal))
			rtd_nam = lvm0['Channel names'][15:20]  # tc14-18
			c0, c1, c2 = lvm0['data'][:, 15], lvm0['data'][:, 16], lvm0['data'][:, 17]
			c3, c4 = lvm0['data'][:, 18], lvm0['data'][:, 19]
			lvmtime, c0, c1, c2, c3, c4 = lvmtime[ii], c0[ii], c1[ii], c2[ii], c3[ii], c4[ii]
			# for validating correct pulse alignment
			# plt.plot(lvmtime, c0)  # ---------------
			# plt.show()  # --------------------------
			ipeak = min(
				[np.where(c0 == max(c0))[0][0], np.where(c1 == max(c1))[0][0], np.where(c2 == max(c2))[0][0],
				 np.where(c3 == max(c3))[0][0], np.where(c4 == max(c4))[0][0]])
			b0, b1, b2, b3, b4 = mean(c0[:ipeak - 5]), mean(c1[:ipeak - 5]), mean(c2[:ipeak - 5]), mean(
				c3[:ipeak - 5]), mean(c4[:ipeak - 5])
			c0, c1, c2, c3, c4 = c0 - b0, c1 - b1, c2 - b2, c3 - b3, c4 - b4
			lvmtime -= lvmtime[ipeak]
			return lvmtime, rtd_nam, c0, c1, c2, c3, c4
	print(f'shot {shot} MISSING from lvm data files')
	return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan


def calorimeter_current_scan_29mar22():
	'''
		plot fwhm vs instantaneous equilibrated dT vs beam current at different voltages
	'''
	
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2)  # 2, figsize=(10, 5))
	
	# direc = 'Y:/thermocouple/Calorimeter/'
	# lvm_files1 = os.listdir(direc)
	# lvm_files1 = [f'{direc}{lvmf}' for lvmf in lvm_files1 if lvmf.startswith('04062022')]  # lvms with date format
	# calibs1 = np.zeros_like(lvm_files1)
	lvm_files, calibs = special_calibration(redo=False)  # figure out calibs per lvm
	# lvm_files = np.append(lvm_files1, lvm_files2)
	# calibs = np.append(calibs1, calibs2)
	
	lbls = ['10', '12.5', '15']
	# no TC data for GUI shots 25, 26, 27, 29, 32, 33, 34, 35
	gui15kv = 508020 + np.array(
		[22, 24, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 49, 51,
		 52])  # NBI GUI shots for 15kv GUI setting plus offset to get to NBI tree shots #s
	star15kv = 508020 + np.array([40, 43, 44, 46, 49])
	gui12p5kv = np.append(508076 + np.array([16, 19, 20, 24, 25, 29, 33, 35, 37, 39]),
	                      508149 + np.array([49, 54, 57]))
	star12p5kv = np.append(508076 + np.array([16, 19, 20, 24, 25, 29, 35, 37]),
	                       508149 + np.array([49, 54, 57]))
	gui10kv = np.append(508149 + np.array([28, 30, 32, 34, 35, 37, 39, 63]),
	                    508268 + np.array([1, 3]))
	star10kv = np.append(508149 + np.array([30, 32, 34, 35, 37, 39, 63]),
	                     508268 + np.array([1, 3]))
	
	# direc = 'Y:/thermocouple/Calorimeter/'
	# lvm_files = os.listdir(direc)
	# lvm_files = [f'{direc}{lvmf}' for lvmf in lvm_files if '2022' in lvmf]  # lvms with date format
	
	stuff = zip([star10kv, star12p5kv, star15kv], [gui10kv, gui12p5kv, gui15kv], lbls)
	for iset, (star, shotset, lbl) in enumerate(stuff):
		fwhm_arr, perv_arr, jtot_arr, ib_arr, dt_arr = np.array([]), np.array([]), np.array([]), np.array([]), np.array(
			[])
		for shot in shotset:
			print(f'analyzing shot {shot}')
			perv_av, jtot, ib_av = avg_perv(shot, Ej=True, ib_return=True)
			perv_arr = np.append(perv_arr, perv_av)
			ib_arr = np.append(ib_arr, ib_av)
			jtot_arr = np.append(jtot_arr, jtot)
			tt, rtd_nam, c0, c1, c2, c3, c4 = special_extraction(lvm_files, shot, calibs)
			if tt is nan:
				fwhm_arr = np.append(fwhm_arr, np.nan)
				dt_arr = np.append(dt_arr, np.nan)
			else:
				dt0, dt1, dt2, dt3, dt4 = max(c0), max(c1), max(c2), max(c3), max(c4)  # already removed offset
				# g = 1/sig*sqrt(2pi)*exp(-x^2/2sig^2)
				r_therm = 45.  # thermocouple distance to center of calorimeter [mm]
				r_edge = 75.  # radius of calorimeter [mm]
				if abs(dt3 - np.mean(
						[dt1, dt2, dt4])) > .25:  # dt3 (upper) is too high, probs insufficiently neutralized- omit
					print(f'insufficient neutralization on shot {shot}')
					sig_arr = [np.nan]
					# sig_arr = [r_therm * np.sqrt(1 / (2 * np.log(dt0 / dt_edge))) for dt_edge in
					#            [dt1, dt2, dt4]]  # [mm]
					dt_arr = np.append(dt_arr, np.nan)
				else:
					sig_arr = [r_therm * np.sqrt(1 / (2 * np.log(dt0 / dt_edge))) for dt_edge in
					           [dt1, dt2, dt3, dt4]]  # [mm]
					islow_decay = np.where((tt > 25) & (tt < 100))[0]
					if len(islow_decay) > 50:
						line2fit = (c0[islow_decay] + c1[islow_decay] + c2[islow_decay] + c3[islow_decay] + c4[
							islow_decay]) / 5.
						time2fit = tt[islow_decay]
						fit = np.polyfit(time2fit, line2fit, 1)
						temp_interp = np.interp(0, tt, tt * fit[0] + fit[1])
						dt_arr = np.append(dt_arr, temp_interp)
					else:
						dt_arr = np.append(dt_arr, np.nan)
				sig = np.mean(sig_arr)
				fwhm = 2 * np.sqrt(2 * np.log(2)) * sig
				power_frac_to_cal = gauss2d_integral(nsig=r_edge / sig)
				fwhm_arr = np.append(fwhm_arr, fwhm)
		
		istar = [np.where(shotset == starr)[0][0] for starr in star]
		ax1.plot(ib_arr, perv_arr * 1.e6, 'o', label=lbl)
		ax1.plot(ib_arr[istar], perv_arr[istar] * 1.e6, 'k+')
		ax1.plot(ib_arr[-1], perv_arr[-1] * 1.e6, 'rx')
		ax2.plot(ib_arr, fwhm_arr, 'o')
		ax2.plot(ib_arr[istar], fwhm_arr[istar], 'k+')
		ax3.plot(ib_arr, dt_arr / jtot_arr, 'o')
		ax3.plot(ib_arr[istar], dt_arr[istar] / jtot_arr[istar], 'k+')
		ax4.plot(ib_arr, dt_arr, 'o')
		ax4.plot(ib_arr[istar], dt_arr[istar], 'k+')
	# for i, sh in enumerate(shotset):
	# 	ax4.annotate(str(sh), (ib_arr[i], dt_arr[i]))
	
	fs = 14
	ax1.legend()
	ax1.set_ylabel('perveance (e-6)', fontsize=fs)
	ax2.set_ylabel('FWHM (mm)', fontsize=fs)
	ax3.set_ylabel('$\Delta$T/$J_{tot}$', fontsize=fs)
	ax4.set_ylabel('cal $\Delta$T', fontsize=fs)
	
	for ax in [ax1, ax2, ax3, ax4]:
		ax.set_xlabel('$I_b$ (A)', fontsize=fs)
		ax.tick_params(labelsize=12)
	plt.tight_layout()


def special_cal_guass_fit_offcenter(shot):
	# center, down, right, up, left
	lvm_files, calibs = special_calibration(redo=False)  # figure out calibs per lvm
	t, rtd_nam, c0, c1, c2, c3, c4 = special_extraction(lvm_files, shot, calibs)
	# t, c0, c1, c2, c3, c4 = cal_temp_sigs(shot, calonly=True)
	dt0, dt1, dt2, dt3, dt4 = max(c0), max(c1), max(c2), max(c3), max(c4)  # already removed offset
	'''
	thermocouples are 45 mm off center
	calorimeter diameter 150 mm
	'''
	maxdt = max([dt0, dt1, dt2, dt3, dt4])
	normvals = [dt0 / maxdt, dt1 / maxdt, dt2 / maxdt, dt3 / maxdt, dt4 / maxdt]  # normalized thermocouple data
	xtherms, ytherms = [0, 0, 45, 0, -45], [0, -45, 0, 45, 0]  # x, y coords of thermocouples [mm]
	
	def minfunc(xysig):
		rtherms = [np.sqrt((xtherms[i] - xysig[0]) ** 2 + (ytherms[i] - xysig[1]) ** 2) for i in range(5)]  # [mm]
		thermvals = [np.exp(-rtherms[i] ** 2 / (2 * xysig[2] ** 2)) for i in range(5)]
		normsim = [thermvals[i] / max(thermvals) for i in range(5)]  # normalize to compare to normalized real data
		diff = np.array([normvals[i] - normsim[i] for i in range(5)]) * 10000
		return np.sqrt(np.sum(diff ** 2))
	
	xguess = [0, 0, 10]  # guess xoffset=0, yoffset=0, sigma=10mm
	res = minimize(minfunc, xguess, method='SLSQP')
	print(res.message)
	
	fit = res.x  # xysig?
	rtherms = [np.sqrt((xtherms[i] - fit[0]) ** 2 + (ytherms[i] - fit[1]) ** 2) for i in range(5)]  # [mm]
	thermvals = [np.exp(-rtherms[i] ** 2 / (2 * fit[2] ** 2)) for i in range(5)]
	normsim = [thermvals[i] / max(thermvals) for i in range(5)]  # normalize to compare to normalized real data
	fwhm = 2 * np.sqrt(2 * np.log(2)) * fit[2]
	
	fig, (ax1, ax2) = plt.subplots(ncols=2)
	for c, l in zip([c0, c1, c2, c3, c4], ['center', 'lower', 'right', 'upper', 'left']):
		ax1.plot(t, c, label=l)
	ax1.legend()
	ax1.set_xlabel('time')
	ax1.set_ylabel('Temp \degree C')
	rr = np.linspace(0, 1.1 * max(rtherms), endpoint=True)
	ax2.plot(rr, np.exp(-rr ** 2 / (2 * fit[2] ** 2)))
	ax2.plot(rtherms, normvals, 's')
	print(f'x-shift: {fit[0]:.1f} mm\ny-shift: {fit[1]:.1f} mm\nfwhm: {fwhm:.1f} mm')
	plt.show()


def cal_guass_fit_offcenter(shot, ret='fwhm', doplot=False):
	# center, down, right, up, left
	t, c0, c1, c2, c3, c4 = cal_temp_sigs(shot, calonly=True)
	dt0, dt1, dt2, dt3, dt4 = max(c0), max(c1), max(c2), max(c3), max(c4)  # already removed offset
	'''
	thermocouples are 45 mm off center
	calorimeter diameter 150 mm
	'''
	maxdt = max([dt0, dt1, dt2, dt3, dt4])
	normvals = [dt0 / maxdt, dt1 / maxdt, dt2 / maxdt, dt3 / maxdt, dt4 / maxdt]  # normalized thermocouple data
	xtherms, ytherms = [0, 0, 45, 0, -45], [0, -45, 0, 45, 0]  # x, y coords of thermocouples [mm]
	
	def minfunc(xysig):
		rtherms = [np.sqrt((xtherms[i] - xysig[0]) ** 2 + (ytherms[i] - xysig[1]) ** 2) for i in range(5)]  # [mm]
		thermvals = [np.exp(-rtherms[i] ** 2 / (2 * xysig[2] ** 2)) for i in range(5)]
		normsim = [thermvals[i] / max(thermvals) for i in range(5)]  # normalize to compare to normalized real data
		diff = np.array([normvals[i] - normsim[i] for i in range(5)]) * 10000
		return np.sqrt(np.sum(diff ** 2))
	
	xguess = [0, 0, 10]  # guess xoffset=0, yoffset=0, sigma=10mm
	res = minimize(minfunc, xguess, method='SLSQP')
	print(res.message)
	
	fit = res.x
	rtherms = [np.sqrt((xtherms[i] - fit[0]) ** 2 + (ytherms[i] - fit[1]) ** 2) for i in range(5)]  # [mm]
	thermvals = [np.exp(-rtherms[i] ** 2 / (2 * fit[2] ** 2)) for i in range(5)]
	normsim = [thermvals[i] / max(thermvals) for i in range(5)]  # normalize to compare to normalized real data
	fwhm = 2 * np.sqrt(2 * np.log(2)) * abs(fit[2])
	updown = dt3 / dt1  # ratio of upper/lower, some measure of neutralization efficiency- perm magnet bends res ion frac up
	if doplot:
		fig, (ax1, ax2) = plt.subplots(ncols=2)
		for c, l in zip([c0, c1, c2, c3, c4], ['center', 'lower', 'right', 'upper', 'left']):
			ax1.plot(t, c, label=l)
		ax1.legend()
		ax1.set_xlabel('time')
		ax1.set_ylabel('Temp \degree C')
		rr = np.linspace(0, 1.1 * max(rtherms), endpoint=True)
		ax2.plot(rr, np.exp(-rr ** 2 / (2 * fit[2] ** 2)))
		ax2.plot(rtherms, normvals, 's')
		print(f'x-shift: {fit[0]:.1f} mm\ny-shift: {fit[1]:.1f} mm\nfwhm: {fwhm:.1f} mm')
		plt.show()
	
	if ret.lower() == 'fwhm':
		return fwhm, updown, fit[0], fit[1]
	# elif ret.lower() == 'power_frac_to_cal':
	# 	return power_frac_to_cal, sufficient_neut
	elif ret.lower() == 'efold':
		return np.sqrt(2) * fit[2], updown, fit[0], fit[1]


def cal_gauss_fit(shot, ret='fwhm'):
	"""
	:param shot:
	:param ret:  'fwhm'- returns full width at half max of guassian
				'power_frac_to_cal'- returns power fraction that hits calorimeter
				'efold'- returns e-folding length of gaussian
	:return:
	"""
	# center, down, right, up, left
	t, c0, c1, c2, c3, c4 = cal_temp_sigs(shot, calonly=True)
	dt0, dt1, dt2, dt3, dt4 = max(c0), max(c1), max(c2), max(c3), max(c4)  # already removed offset
	'''
	thermocouples are 45 mm off center
	calorimeter diameter 150 mm
	'''
	# g = 1/sig*sqrt(2pi)*exp(-x^2/2sig^2)
	r_therm = 45.  # thermocouple distance to center of calorimeter [mm]
	r_edge = 75.  # radius of calorimeter [mm]
	if abs(dt3 - np.mean([dt1, dt2, dt4])) > .25:  # dt3 (upper) is too high, probably insufficiently neutralized- omit
		print(f'insufficient neutralization for shot {shot}')
		sufficient_neut = 0
		sig_arr = [r_therm * np.sqrt(1 / (2 * np.log(dt0 / dt_edge))) for dt_edge in [dt1, dt2, dt4]]  # [mm]
	else:
		sufficient_neut = 1
		sig_arr = [r_therm * np.sqrt(1 / (2 * np.log(dt0 / dt_edge))) for dt_edge in [dt1, dt2, dt3, dt4]]  # [mm]
	sig = np.mean(sig_arr)
	fwhm = 2 * np.sqrt(2 * np.log(2)) * sig
	power_frac_to_cal = gauss2d_integral(nsig=r_edge / sig)
	if ret.lower() == 'fwhm':
		return fwhm, sufficient_neut
	elif ret.lower() == 'power_frac_to_cal':
		return power_frac_to_cal, sufficient_neut
	elif ret.lower() == 'efold':
		return np.sqrt(2) * sig, sufficient_neut


def gauss2d_integral(nsig=1., inspect=False):
	# nsig: compute volume contained within nsig of 2d gaussian
	x = np.linspace(-100, 100, num=1000, endpoint=True)
	dx, sig = x[1] - x[0], 10
	y = 1 / (sig * np.sqrt(2 * np.pi)) * np.exp(-x ** 2 / (2 * sig ** 2))
	area = np.sum(y) * dx  # check: this should be = 1
	vol2d = np.sum(y * 2 * np.pi * abs(x) * dx) / 2
	vol2d_calc = np.sqrt(2 * np.pi) * sig
	y2 = np.copy(y)
	y2[abs(x) > nsig * sig] = 0
	inside_vol = np.sum(y2 * 2 * np.pi * abs(x) * dx) / 2  # double counting since x extends both pos and neg
	ratio_inside = inside_vol / vol2d
	if inspect:
		plt.plot(x, y)
		plt.show()
	return ratio_inside


def neutralization_scan_7feb22():
	desc = ['10', '20', '30']
	sh20psi = 105000 + np.array([81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100])
	# starts at 18 were actually neutralizer=off, since beam lasts until 15, should be same as starting at 16
	ns20psi = np.array([18, 18, 16, 16, 14, 14, 12, 12, 10, 10, 8, 8, 6, 6, 4, 4, 2, 2, 0, 0])
	sh30psi = 105100 + np.array([1, 2, 3, 4, 5, 7, 8, 9, 10])
	ns30psi = np.array([8, 8, 6, 10, 4, 12, 2, 2, 14])
	sh10psi = 105100 + np.array([11, 12, 13, 14, 15, 16, 17, 18])
	ns10psi = np.array([8, 6, 10, 4, 12, 2, 14, 14])
	fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5))
	for ii, (ns, ss) in enumerate(zip([ns10psi, ns20psi, ns30psi], [sh10psi, sh20psi, sh30psi])):
		dt_meas, dt_pred, pow_frac, dt_up = np.array([]), np.array([]), np.array([]), np.array([])
		for shot in ss:
			dtm, dtp, dt, _ = cal_dtemp(shot, dbug=False, more=True)
			pf, _ = cal_gauss_fit(shot)
			dt_meas = np.append(dt_meas, dtm)
			dt_pred = np.append(dt_pred, dtp)
			pow_frac = np.append(pow_frac, pf)
			dt_up = np.append(dt_up, dt[3] - np.mean([dt[1], dt[2], dt[4]]))  # up - mean(down, right, left)
		ax1.plot(ns, dt_meas / (pow_frac * dt_pred), 'o', label=desc[ii])
		ax2.plot(ns, dt_up, 'o')
	ax1.legend(title='GVC psi')
	xt = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]
	for ax in [ax1, ax2]:
		ax.set_xticks(xt)
		ax.set_xlabel('neutralizer start (ms)')
	ax1.set_ylabel('meas/pred power fraction')
	ax2.set_ylabel('extra $\Delta T$ on upper')
	plt.tight_layout()


def perveance_scan_28feb22():
	'''
	HAL data today: compare to perveance_scan_7feb22()
	'''
	sh_v = 105300 + np.array([40, 41, 42, 44, 45, 53, 54, 55, 57, 58, 60, 61])  # vb scan
	sh_i = 105300 + np.array([65, 66, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 83, 84, 85, 87])  # ib scan
	sh3 = 105300 + np.array([89, 90, 93, 94])  # higher current, pulling down voltage to access higher perv
	sh4 = 105300 + np.array(
		[95, 96, 97, 98, 99, 101, 102, 103, 104, 105, 106, 107, 108])  # Iarc=400, voltage scan at lower ib
	fig, (ax1) = plt.subplots(ncols=1)  # 2, figsize=(10, 5))
	for i, shotset in enumerate([sh_v, sh_i, sh3, sh4]):
		print(shotset)
		perv_arr = np.array([])
		for shot in shotset:
			perv_arr = np.append(perv_arr, avg_perv(shot))
		yvals = np.ones_like(perv_arr) + i / 10.
		ax1.plot(perv_arr * 1.e6, yvals, 'o')
		ax1.plot(perv_arr[-1] * 1.e6, yvals[-1], 'kd')
	
	ax1.set_xlabel('perveance (e-6)')
	ax1.set_ylabel('placeholder')
	plt.tight_layout()


def stray_field_test_6apr22():
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	tp = np.linspace(-50, 100, endpoint=True)
	fig = plt.figure(tight_layout=True, figsize=(10, 8))
	gs = gridspec.GridSpec(2, 3)
	ax1 = fig.add_subplot(gs[0, 0])
	ax2 = fig.add_subplot(gs[0, 1])
	ax3 = fig.add_subplot(gs[0, 2])
	ax_bot = fig.add_subplot(gs[1, :])
	topax = [ax1, ax2, ax3]
	
	lbls = ['no fields', 'high fields', 'highest fields']
	nofields = 105500 + np.array([94, 95, 96])
	fields = 105500 + np.array([88, 89, 93])
	highfields = 105500 + np.array([97, 99, 100])
	dt_arr = np.zeros((5, 3))
	err_arr = np.zeros((5, 3))
	for i, (lbl, shots) in enumerate(zip(lbls, [nofields, fields, highfields])):
		etot = 0.
		for sh in shots:
			_, shot_etot = avg_perv(sh, Ej=True)
			etot += shot_etot
		etot = etot / len(shots) * 1.e-3  # kJ
		cntr, dwn, rght, uppr, lft, shot_dt, minmax = get_thermocouple_numbers(shots, tp, ret_minmax=True)
		errs = (minmax[:, 1] - minmax[:, 0]) / 2.
		err_arr[:, i] = errs
		dt_norm = [dt / etot for dt in shot_dt]
		dt_arr[:, i] = dt_norm  # shot_dt
		topax[i].plot(tp, cntr / etot, clrs[0], label='center')  # correct labeling
		topax[i].plot(tp, dwn / etot, clrs[1], label='down')
		topax[i].plot(tp, rght / etot, clrs[2], label='right')
		topax[i].plot(tp, uppr / etot, clrs[3], label='up')
		topax[i].plot(tp, lft / etot, clrs[4], label='left')
	
	for i in range(5):
		ax_bot.errorbar(np.array([-1, 0, 1]) - .1 + .05 * i, dt_arr[i, :], yerr=err_arr[i, :], c=clrs[i], capsize=5)
	
	fs = 12
	ax1.legend(fontsize=fs)
	ax1.set_ylabel('thermocouple\n$\Delta T/E_{tot}$ (\degree C/kJ)', fontsize=fs)
	# ax1.set_xlabel('time rel to beam')
	ax2.set_ylim(ax1.get_ylim())
	ax2.set_xlabel('time rel to beam', fontsize=fs)
	ax_bot.set_ylabel('thermocouple\n$\Delta T/E_{tot}$ (deg C/kJ)', fontsize=fs)
	ax_bot.set_xlim((-1.5, 1.5))
	ax_bot.set_xticks([-1, 0, 1])
	
	ax1.set_title('No Fields', fontsize=fs)
	ax2.set_title('With Fields', fontsize=fs)
	ax3.set_title('High Fields', fontsize=fs)
	ax_bot.set_xticklabels(['no fields', 'fields', 'high fields'], fontsize=fs)
	for ax in [ax1, ax2, ax3, ax_bot]:
		ax.tick_params(labelsize=fs)
	plt.tight_layout()
	plt.show()


def perveance_scan_29Jun22(plot_vs='fwhm'):
	"""
	beam misaligned on calorimeter, rtan~25cm
	looking to see if beam is still neutralizing and getting good power to calorimeter
	"""
	if plot_vs.lower() == 'fwhm' or plot_vs.lower() == 'efold':
		units = ' (mm)'
	else:
		units = ''
	
	fig, (ax1) = plt.subplots(ncols=1)  # 2, figsize=(10, 5))
	sh1 = 106600 + np.array([60, 62, 64, 65, 67, 68, 69, 70, 71, 72])  # Current scan at 15kV, rtan~25cm
	sh2 = 106600 + np.array([74, 76, 78, 79])  # Current scan at 12.5kV, rtan~25cm
	sh3 = 106700 + np.array([4, 8, 10, 11, 12, 13])  # Current scan at 15kV, rtan~21cm, source centered on neutralizer
	sh4 = 106700 + np.array(
		[17, 20, 23, 25, 29, 30])  # Current scan at 12.5kV, rtan~21cm, source centered on neutralizer
	descs = ['15kV, rtan=25cm', '12.5kV, rtan=25cm', '15kV, rtan=21cm', '12.5kV, rtan=21cm']
	for sh in [sh1, sh2, sh3, sh4]:
		fwhm_arr, perv_arr = np.array([]), np.array([])
		updown_arr, xshift_arr, yshift_arr = np.array([]), np.array([]), np.array([])
		for ish, shot in enumerate(sh):
			perv_arr = np.append(perv_arr, avg_perv(shot))
			fwhm, updown, xshift, yshift = cal_guass_fit_offcenter(shot, ret=plot_vs, doplot=False)
			fwhm_arr = np.append(fwhm_arr, fwhm)
			updown_arr = np.append(updown_arr, updown)
			xshift_arr = np.append(xshift_arr, xshift)
			yshift_arr = np.append(yshift_arr, yshift)
		fs = 14
		ax1.plot(perv_arr * 1.e6, fwhm_arr, 'o')
		for i in range(len(xshift_arr)):
			note = f'{xshift_arr[i]:.1f}\n{yshift_arr[i]:.1f}'
			ax1.annotate(note, (perv_arr[i] * 1.e6, fwhm_arr[i]))
	ax1.set_xlabel('perveance ($\\times 10^{-6}$)', fontsize=fs)
	ax1.set_ylabel(f'{plot_vs}{units}', fontsize=fs)
	ax1.tick_params(labelsize=fs)
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	for i, desc in enumerate(descs):
		ax1.plot(np.nan, np.nan, 'o', c=clrs[i], label=desc)
	ax1.legend(loc='upper right')
	plt.tight_layout()


def power_to_calorimeter_check_26aug22():
	shots13kv = 107000 + np.array([294, 296, 301, 302, 303, 304, 315, 316, 317])
	ms13kv = np.array([5, 5, 5, 7.5, 7.5, 7.5, 10, 10, 10])
	shots8kv = 107300 + np.array([5, 7, 9, 12, 13, 14, 18, 19, 21])
	ms8kv = np.array([7.5, 7.5, 7.5, 10, 10, 10, 5, 5, 5])
	desc = ['13 kV', '8 kV']
	fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5), sharey='row')
	for ii, (ss, pl) in enumerate(zip([shots13kv, shots8kv], [ms13kv, ms8kv])):
		dt_meas, dt_pred = np.array([]), np.array([])
		fwhm_arr, xshift_arr, yshift_arr = np.array([]), np.array([]), np.array([])
		for shot in ss:
			# print(f'{shot}')
			dtm, dtp, dt, _ = cal_dtemp(shot, dbug=False, more=True)
			fwhm, updown, xshift, yshift = cal_guass_fit_offcenter(shot, ret='fwhm', doplot=False)
			fwhm_arr = np.append(fwhm_arr, fwhm)
			xshift_arr = np.append(xshift_arr, xshift)
			yshift_arr = np.append(yshift_arr, yshift)
			dt_meas = np.append(dt_meas, dtm)
			dt_pred = np.append(dt_pred, dtp)
		power_frac = dt_meas / dt_pred
		ax1.plot(pl, power_frac, 'o', label=desc[ii])
		ax2.plot(fwhm_arr, power_frac, 'o', label=desc[ii])
		for i in range(len(xshift_arr)):
			note = f'{xshift_arr[i]:.1f}\n{yshift_arr[i]:.1f}'
			ax1.annotate(note, (pl[i], power_frac[i]))
	ax1.legend()
	ax1.set_ylabel('meas/pred power fraction')
	ax1.set_xlabel('pulse length (ms)')
	ax2.set_xlabel('fwhm (mm)')
	ax1.set_xlim((4, 11))
	plt.tight_layout()


if __name__ == '__main__':
	# perveance_scan_7feb22(plot_vs='FWHM')
	# new_valve_data()
	# new_valve_cathode()
	# new_valve_neutralizer()
	# offcenter = 508397
	# centered = 508188
	# cal_guass_fit_offcenter(offcenter)
	# calorimeter_current_scan_29mar22()
	# calorimeter_21Jan()  # max 53%
	# stray_field_test_6apr22()
	# calorimeter_current_scan_29mar22()
	
	# compare_neutralizer_onoff()
	# perveance_scan_7feb22(plot_vs='fwhm')
	# perveance_scan_28feb22()
	# neutralization_scan_7feb22()
	# realigned_beam_4feb22()
	# gauss2d_integral(2, inspect=True)
	# cal_gauss_fit([105048])  # 65
	# calorimeter_positional_variation()
	beam_realignment()
	# calorimeter_21Jan()  # coincides w/settings used for 14Jan Perveance scan
	# perveance_scan_29Jun22()
	# power_to_calorimeter_check_26aug22()
	plt.show()
