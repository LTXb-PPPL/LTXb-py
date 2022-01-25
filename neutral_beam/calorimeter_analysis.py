from neutral_beam.rtd_analysis import avg_perv
import matplotlib.pyplot as plt
import numpy as np
from helpful_stuff import SimpleSignal
import lvm_read
from bills_LTX_MDSplus_toolbox import *
import os

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
		perv, pwr = ones_like(vb), ones_like(vb)
		perv[:], pwr[:] = np.nan, np.nan
		perv[t_window] = ib[t_window] / vb[t_window] ** 1.5
		pwr[t_window] = ib[t_window] * vb[t_window]  # [W]
		av_perv, av_pwr = np.nanmean(perv), np.nanmean(pwr)
		dperv, dpwr = (np.nanmax(perv) - np.nanmin(perv)) / 2., (np.nanmax(pwr) - np.nanmin(pwr)) / 2.
		return av_perv, dperv, av_pwr, dpwr
	except mds.mdsExceptions.TreeNODATA:
		print(f'trouble fetching MDSPlus data for shot {shot}')
		return np.nan, np.nan, np.nan, np.nan


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
		lvm0 = lvm[0]
		time = lvm0['data'][:, 0]
		rtd_nam = lvm0['Channel names'][4:9]  # tc14-18
		
		c0, c1, c2 = lvm0['data'][:, 4], lvm0['data'][:, 5], lvm0['data'][:, 6]
		c3, c4 = lvm0['data'][:, 7], lvm0['data'][:, 8]
		it = np.where(time < 300)
		time, c0, c1, c2, c3, c4 = time[it], c0[it], c1[it], c2[it], c3[it], c4[it]
		ipeak = np.where(c0 == max(c0))[0][0]
		baseline = mean((c0[:ipeak - 1] + c1[:ipeak - 1] + c2[:ipeak - 1] + c3[:ipeak - 1] + c4[:ipeak - 1]) / 5.)
		islow_decay = np.where((time > time[ipeak] + 25) & (time < time[ipeak] + 100))[0]
		if len(islow_decay) < 100:
			print(f'--insufficient data to analyze for shot {shot} in file {lvm_file}')
			return nodata
		else:
			dt0, dt1, dt2, dt3, dt4 = c0[ipeak] - mean(c0[0:ipeak - 5]), c1[ipeak] - mean(c1[0:ipeak - 5]), c2[
				ipeak] - mean(c2[0:ipeak - 5]), c3[ipeak] - mean(c3[0:ipeak - 5]), c4[ipeak] - mean(c4[0:ipeak - 5])
			line2fit = (c0[islow_decay] + c1[islow_decay] + c2[islow_decay] + c3[islow_decay] + c4[islow_decay]) / 5.
			time2fit = time[islow_decay]
			fit = np.polyfit(time2fit, line2fit, 1)
			temp_interp = np.interp(time[ipeak], time, time * fit[0] + fit[1])
			meas = temp_interp - baseline
			pred = tot_joules / c_ps / m_copp  # predicted dT (degC) from injected beam energy
			
			if dbug:
				plt.suptitle(f'{shot}')
				plt.plot(time, c0, label='tc0')
				# plt.axvline(time[islow_decay[0]])
				# plt.axvline(time[islow_decay[-1]])
				plt.plot(time, c1, label='tc1')
				plt.plot(time, c2, label='tc2')
				plt.plot(time, c3, label='tc3')
				plt.plot(time, c4, label='tc4')
				plt.legend()
				plt.xlim((0, 60))
				# plt.axhline(baseline)
				# plt.plot(time, time * fit[0] + fit[1])
				plt.xlabel('time (s)')
				plt.ylabel('temp ($\degree$C)')
				plt.show()
			
			return meas, pred, [dt0, dt1, dt2, dt3, dt4], [ipfit[0], stdev]


def cal_temp_sigs(shot):
	nodata = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
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
		lvm0 = lvm[0]
		time = lvm0['data'][:, 0]
		rtd_nam = lvm0['Channel names'][4:9]  # tc14-18
		
		c0, c1, c2 = lvm0['data'][:, 4], lvm0['data'][:, 5], lvm0['data'][:, 6]
		c3, c4 = lvm0['data'][:, 7], lvm0['data'][:, 8]
		it = np.where(time < 300)
		time, c0, c1, c2, c3, c4 = time[it], c0[it], c1[it], c2[it], c3[it], c4[it]
		ipeak = np.where(c0 == max(c0))[0][0]
		b0, b1, b2, b3, b4 = mean(c0[:ipeak - 1]), mean(c1[:ipeak - 1]), mean(c2[:ipeak - 1]), mean(
			c3[:ipeak - 1]), mean(c4[:ipeak - 1])
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


def compare_thermocouple_signals(beforenafter=True):
	import matplotlib.gridspec as gridspec
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	
	# new orifice, original alignment
	shots21jan22 = 104800 + np.array([66, 68, 69, 70, 71, 72, 74, 75, 76, 77, 78, 80, 85, 86, 87, 89, 92, 93])
	# new orifice, new alignment
	shots25jan22 = 104900 + np.array([26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 43, 44, 45, 46])
	# new orifice, new alignment, no neutralizer
	shots25jan22_noneut = 104900 + np.array([38, 41, 48, 49, 50, 51])
	tp = np.linspace(-50, 100, endpoint=True)
	
	def get_numbers(shots):
		c0p, c1p, c2p, c3p, c4p = np.zeros_like(tp), np.zeros_like(tp), np.zeros_like(tp), np.zeros_like(
			tp), np.zeros_like(tp)
		num = 0.
		for shot in shots:
			if 12.5 < avg_perv(shot) * 1.e6 < 18:
				num += 1
				t, c0, c1, c2, c3, c4 = cal_temp_sigs(shot)
				c0p += np.interp(tp, t, c0)
				c1p += np.interp(tp, t, c1)
				c2p += np.interp(tp, t, c2)
				c3p += np.interp(tp, t, c3)
				c4p += np.interp(tp, t, c4)
		c0p, c1p, c2p, c3p, c4p = c0p / num, c1p / num, c2p / num, c3p / num, c4p / num
		ipeak = np.where(c0p == max(c0p))[0][0]
		# dt0, dt1, dt2, dt3, dt4 = max(c0p), max(c1p), max(c2p), max(c3p), max(c4p)
		dt0, dt1, dt2, dt3, dt4 = c0p[ipeak], c1p[ipeak], c2p[ipeak], c3p[ipeak], c4p[ipeak]
		return c0p, c1p, c2p, c3p, c4p, [dt0, dt1, dt2, dt3, dt4]
	
	if beforenafter:
		pre0, pre1, pre2, pre3, pre4, dt_pre = get_numbers(shots21jan22)
		post0, post1, post2, post3, post4, dt_post = get_numbers(shots25jan22)
	else: # neut vs no neutralizer
		pre0, pre1, pre2, pre3, pre4, dt_pre = get_numbers(shots25jan22)
		post0, post1, post2, post3, post4, dt_post = get_numbers(shots25jan22_noneut)

	fig = plt.figure(tight_layout=True, figsize=(10, 8))
	gs = gridspec.GridSpec(2, 2)
	ax1 = fig.add_subplot(gs[0, 0])
	ax2 = fig.add_subplot(gs[0, 1])
	ax3 = fig.add_subplot(gs[1, :])
	ax1.plot(tp, pre0, clrs[0], label='center')
	ax1.plot(tp, pre1, clrs[1], label='right (down)')
	ax1.plot(tp, pre2, clrs[2], label='up (right)')
	ax1.plot(tp, pre3, clrs[3], label='left (up)')
	ax1.plot(tp, pre4, clrs[4], label='down (left)')
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

	if beforenafter:
		ax1.set_title('Before')
		ax2.set_title('After')
		ax3.set_xlabel('beam realignment')
		ax3.set_xticklabels(['before', 'after'])
	else:
		ax1.set_title('Neutralizer ON')
		ax2.set_title('Neutralizer OFF')
		ax3.set_xticklabels(['neut ON', 'neut OFF'])

	plt.show()


# a = 1
#
# offs = [0, 10, 0]
# # [old orifices & aiming, new orifices & old aiming, new orifices & aiming]
# for i, shot in enumerate([506587, 104866, 104936]):
# 	t, c0, c1, c2, c3, c4 = cal_temp_sigs(shot)
# 	plt.plot(t + offs[i], c0, clrs[0])
# 	plt.plot(t + offs[i], c1, clrs[1])
# 	plt.plot(t + offs[i], c2, clrs[2])
# 	plt.plot(t + offs[i], c3, clrs[3])
# 	plt.plot(t + offs[i], c4, clrs[4])
# plt.xlim((-5, 50))
# plt.xlabel('time')
# plt.ylabel('temp increase')
# plt.show()


if __name__ == '__main__':
	compare_thermocouple_signals(beforenafter=True)
	# calorimeter_21Jan()  # coincides w/settings used for 14Jan Perveance scan
	plt.show()
# big_scan = 0
# pulse_timing = 0
# anode_cathode_puffing = 0
# new_valve_data = 0
# new_valve_neutralizer = 0
# new_valve_cathode = 0
