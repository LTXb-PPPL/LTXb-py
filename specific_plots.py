import matplotlib.pyplot as plt
import numpy as np
from helpful_stuff import get_tree_conn, get_data, make_patch_spines_invisible, SimpleSignal, is_nbi_shot
from neutral_beam.calorimeter_analysis import cal_temp_sigs, cal_dtemp


def shot_analysis_01mar22():
	shots_1mar22 = np.arange(33) + 105426  # shots 105177-105201
	find_good_nbi = 0
	if find_good_nbi:
		good_nbi = []
		for shot in shots_1mar22:
			t = get_tree_conn(shot, treename='ltx_b')
			if is_nbi_shot(shot, t):
				good_nbi.append(shot)
	else:
		good_nbi = [105426, 105427, 105428, 105429, 105430, 105431, 105432, 105433, 105435, 105436, 105437, 105438,
		            105439, 105440, 105441, 105442, 105443, 105444, 105450, 105451, 105452, 105453, 105454, 105458]
	
	picked_shots = 105400 + np.array([32, 30, 36, 28, 40])
	t460 = [26, 31, 32]  # 32 best
	t463 = [30, 51, 52, 53, 54, 58]  # not 58, 30 is best
	t465 = [29, 35, 36, 41, 42, 43, 44]  # not 41, pick 44 or average
	t468 = [28]
	t470 = [37, 38, 39, 40]  # all good, pick 40 or average
	
	# for shot in t460:
	# 	t = get_tree_conn(shot+105400, treename='ltx_b')
	# 	(tnel, nel) = get_data(t, '.diagnostics.microwave.interferom.phase_comp.ne_l')  # [m^-2]
	# 	plt.plot(tnel, nel, label=shot)
	# t = get_tree_conn(105428, treename='ltx_b')
	# (tnel, nel) = get_data(t, '.diagnostics.microwave.interferom.phase_comp.ne_l')  # [m^-2]
	# plt.plot(tnel, nel, 'k')
	# plt.xlim((.445, .485))
	# plt.legend()
	
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	fig, ax = plt.subplots()
	axr = ax.twinx()
	axrr = ax.twinx()
	axrr.spines["right"].set_position(("axes", 1.2))
	make_patch_spines_invisible(axrr)
	axrr.spines["right"].set_visible(True)
	axrr.yaxis.set_label_position('right')
	axrr.yaxis.set_ticks_position('right')
	for i, shot in enumerate(picked_shots):
		try:
			t = get_tree_conn(shot, treename='ltx_b')
			(tip, ip) = get_data(t, '.diagnostics.magnetic.ip_rog')
			(tnel, nel) = get_data(t, '.diagnostics.microwave.interferom.phase_comp.ne_l')  # [m^-2]
			(tibeam, ibeam) = get_data(t, '.oper_diags.ltx_nbi.source_diags.i_hvps')
			(tbeam, vbeam) = get_data(t, '.oper_diags.ltx_nbi.source_diags.v_hvps')
			ibeam = np.interp(tbeam, tibeam, ibeam)
			pbeam = ibeam * vbeam
			ax.plot(tip * 1e3, -ip * 1.e-3, c=clrs[i], label=shot)
			axr.plot(tnel * 1e3, nel / .5 * 1.e-19, c=clrs[i])  # divide by 0.5m to get <ne>
			axrr.plot(tbeam * 1e3, pbeam * 1.e-6, c=clrs[i])
		except:
			print(f'problem with {shot}')
	ax.set_xlim((445, 485))
	ax.set_ylim((0, 130))
	axr.set_ylim((0, 1.5))
	axrr.set_ylim(0, 2)
	fs = 12
	ax.set_xlabel('time (ms)', fontsize=fs + 2)
	ax.set_ylabel('$I_p$ (kA)', fontsize=fs + 2)
	axr.set_ylabel('<$n_e$> (10$^{19}$ m$^{-3}$)', fontsize=fs + 2)
	axrr.set_ylabel('$P_{nbi}$ (MW)', fontsize=fs + 2)
	for aax in [ax, axr, axrr]:
		aax.tick_params(labelsize=fs)
	ax.legend(fontsize=fs, loc='upper left')
	plt.tight_layout()
	plt.show()


def beam_nobeam_103446_465():
	twin = [450, 474]
	fig, ax = plt.subplots()
	axr = ax.twinx()
	for clr, shot in zip(['b', 'r'], [103465, 103446]):
		try:
			tree = get_tree_conn(shot, treename='ltx_b')
			(tip, ip) = get_data(tree, '.diagnostics.magnetic.ip_rog', times=None)
			(tnel, nel) = get_data(tree, '.diagnostics.microwave.interferom.phase_comp.ne_l')  # [m^-2]
			ax.plot(tip * 1e3, -ip * 1.e-3, label=shot, c=clr)
			axr.plot(tnel * 1e3, nel / .5 * 1.e-19, c=clr)  # divide by 0.5m to get <ne>
			if shot == 103465:
				(tibeam, ibeam) = get_data(tree, '.oper_diags.ltx_nbi.source_diags.i_hvps')
				(tbeam, vbeam) = get_data(tree, '.oper_diags.ltx_nbi.source_diags.v_hvps')
				ibeam = np.interp(tbeam, tibeam, ibeam)
				axr.plot(tbeam * 1e3, ibeam * vbeam * 1.e-6, c=clr)
		except:
			print(f'shot {shot} did not pan out')
	fs = 12
	ax.set_xlim(twin)
	ax.set_ylim((0, 150))
	ax.set_ylabel('$I_p$ (kA)', fontsize=fs)
	ax.set_xlabel('time (ms)', fontsize=fs)
	axr.set_ylabel('<$n_e$> (10$^{19}$ m$^{-3}$) $P_{nbi}$ (MW)', fontsize=fs)
	axr.set_ylim((0, 1.5))
	ax.legend(fontsize=fs)
	ax.tick_params(labelsize=fs)
	axr.tick_params(labelsize=fs)
	plt.tight_layout()
	plt.show()


def compare_thermocouple_sigs_to_russian_trace():
	prealigned = 104866
	bestshot = 105124  # run section below to set bestshot
	# sh1 = 105100 + np.array(
	# 	[24, 27, 31, 32, 35, 45, 47, 48, 49, 50, 51, 52, 53, 54, 55, 57, 58, 59, 61, 62, 63, 64, 65, 67, 68, 69, 71, 73,
	# 	 74, 75, 76])
	# max_dt, bestshot = -1, -1
	# dt_arr = []
	# for shot in sh1:
	# 	t, c0, c1, c2, c3, c4 = cal_temp_sigs(shot, calonly=True)
	# 	dt_arr.append(max(c0))
	# 	if max(c0) > max_dt:
	# 		max_dt = max(c0)
	# 		bestshot = shot
	# print(f'best shot is {bestshot} with dt={max_dt}')
	pret, pre0, pre1, pre2, pre3, pre4 = cal_temp_sigs(prealigned, calonly=True)
	t, c0, c1, c2, c3, c4 = cal_temp_sigs(bestshot, calonly=True)
	
	fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
	lbls = ['center', 'lower', 'right', 'upper', 'left']
	for i, c in enumerate([pre0, pre1, pre2, pre3, pre4]):
		ax1.plot(pret, c, label=lbls[i])
	for c in [c0, c1, c2, c3, c4]:
		ax2.plot(t, c)
	fs = 12
	ax2.set_xlabel('time (s)', fontsize=fs)
	ax1.set_ylabel('temp ($\degree C$)', fontsize=fs)
	ax2.set_ylabel('temp ($\degree C$)', fontsize=fs)
	# ax1.legend(title='thermocouple', fontsize=fs, ncol=2)
	for ax in [ax1, ax2]:
		ax.tick_params(labelsize=fs)
		ax.grid()
	plt.tight_layout()
	plt.xlim(-5, 30)
	plt.show()


def proposal_nfi():
	transp_run = 1000021501  # "optimal" 155kA, 2.3e19 ne(0)
	transp_run = 1057950301  # shot #105795, ~120kA, ~3e19 ne(0)
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	bdens = SimpleSignal(transp_run, '\\bdens')
	bdens2 = np.max(bdens.data, axis=1)  # core nfi (cm^-3)
	ne = SimpleSignal(transp_run, '\\ne')
	ne0 = np.max(ne.data)  # core ne vs time (cm^-3)
	fi_frac = bdens2 / ne0
	fig, ax = plt.subplots()
	axr = ax.twinx()
	ax.plot(bdens.dim2 * 1000., bdens2 * 1.e-12, c=clrs[0])
	axr.plot(bdens.dim2 * 1000., fi_frac, c=clrs[1])
	fs = 12
	ax.set_xlabel('time (ms)', fontsize=fs)
	ax.set_ylabel('$n_{fi}$ 10$^{12}$ cm$^{-3}$', fontsize=fs)
	axr.set_ylabel('$n_{fi}/n_e(0)$', fontsize=fs)
	# ax.spines['left'].set_color(clrs[0])
	# axr.spines['right'].set_color(clrs[1])
	# ax.tick_params(axis='y', colors=clrs[0])
	# axr.tick_params(axis='y', colors=clrs[1])
	for a in [ax, axr]:
		a.tick_params(labelsize=fs)
		a.set_xlim((455, 485))
	axylim = ax.get_ylim()
	ax.fill_between([468, 475.5], axylim[0], axylim[1], facecolor='k', alpha=0.25)
	ax.annotate('beam on', [468, axylim[1] * .95], fontsize=fs)
	ax.set_ylim(axylim)
	# fig.suptitle('discharge #105795', fontsize=fs)
	plt.tight_layout()
	plt.show()


def neutralization_improvements():
	original = 104000 + np.array([55, 56, 58, 60, 67, 68, 69, 71, 74, 75, 76, 77, 78, 79])
	original = np.append(original, np.array(
		[506587, 506589, 506590, 506591, 506592, 506593, 506594, 506596, 506598, 506600, 506601, 506603, 506609]))
	
	searching = np.array([506610, 506612, 506613, 506614, 506616, 506617, 506619, 506620, 506621, 506622, 506626,
	                      506627, 506628, 506630, 506632, 506633, 506634, 506636, 506637, 506638, 506639, 506642,
	                      506643, 506644, 506646, 506647, 506648])
	dt_pred_orig, dt_meas_orig = np.array([]), np.array([])  # shot arrays for joules
	for shot in original:
		dtm, dtp = cal_dtemp(shot, dbug=False)
		dt_meas_orig = np.append(dt_meas_orig, dtm)
		dt_pred_orig = np.append(dt_pred_orig, dtp)
	
	dt_pred_srch, dt_meas_srch = np.array([]), np.array([])  # shot arrays for joules
	for shot in searching:
		dtm, dtp = cal_dtemp(shot, dbug=False)
		dt_meas_srch = np.append(dt_meas_srch, dtm)
		dt_pred_srch = np.append(dt_pred_srch, dtp)
	
	slope_orig = np.nanmean(dt_meas_orig / dt_pred_orig)
	slope_srch = dt_meas_srch[np.where(searching == 506648)[0][0]] / dt_pred_srch[np.where(searching == 506648)[0][0]]
	
	fs = 14
	fig, ax1 = plt.subplots()  # , figsize=(12, 6))
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	ax1.plot([0, max(dt_pred_orig)], [0, max(dt_pred_orig) * slope_orig], '--', c=clrs[0])
	ax1.plot(dt_pred_orig, dt_meas_orig, 'o', c=clrs[0], label='baseline data')
	ax1.annotate(f'{slope_orig * 100.:.2f}%', (max(dt_pred_orig), max(dt_pred_orig) * slope_orig), c=clrs[0],
	             fontsize=fs - 2)
	
	ax1.plot([0, max(dt_pred_srch)], [0, max(dt_pred_srch) * slope_srch], '--', c=clrs[3])
	ax1.plot(dt_pred_srch, dt_meas_srch, 'o', c=clrs[3], label='optimization progress')
	ax1.annotate(f'{slope_srch * 100.:.2f}%', (max(dt_pred_srch), max(dt_pred_srch) * slope_srch), c=clrs[3],
	             fontsize=fs - 2)
	
	# ax1.plot(np.linspace(0, max(dt_pred)), np.linspace(0, max(dt_pred)), 'k--', alpha=.5)
	ax1.set_xlabel('predicted $\Delta T$ ($\degree$C)', fontsize=fs)
	ax1.set_ylabel('measured $\Delta T$ ($\degree$C)', fontsize=fs)
	ax1.tick_params(labelsize=fs)
	ax1.set_xlim(xmax=6.25)
	# ax1.set_title('Calorimeter')
	ax1.legend(fontsize=fs)
	plt.tight_layout()
	plt.show()


if __name__ == '__main__':
	# beam_nobeam_103446_465()
	# shot_analysis_01mar22()
	# compare_thermocouple_sigs_to_russian_trace()
	# proposal_nfi()
	neutralization_improvements()
