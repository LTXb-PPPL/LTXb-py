import matplotlib.pyplot as plt
import numpy as np
from helpful_stuff import get_tree_conn, get_data, make_patch_spines_invisible


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


if __name__ == '__main__':
	# beam_nobeam_103446_465()
	shot_analysis_01mar22()
