"""
To look for more signals, remote onto NoMachine and run traverser
open tree "ltx_b" and shot 101602 (for example)
then navigate around to see which signals are available to plot
"""
from helpful_stuff import smooth, get_current_ltx_shot, is_nbi_shot, is_good_shot
from bills_LTX_MDSplus_toolbox import *
from shotgroup_database import get_shots
import matplotlib.pyplot as plt
import MDSplus


def plot_ip(ax, tree, lbl=None):
	(times, ip) = get_data(tree, '.diagnostics.magnetic.ip_rog', times=None)
	ii = np.where((times >= 0.44) & (times <= 0.48))
	ax.plot(times[ii], -ip[ii], label=lbl)


def plot_sig(ax, tree, node, lbl=None):
	try:
		(times, sig) = get_data(tree, node, times=None)
		ax.plot(times, sig, label=lbl)
	except:
		print('no data for node {}'.format(node))


def avg_density_during_beam(shot):
	try:
		tree = get_tree_conn(shot, treename='ltx_b')
		(tt, ne) = get_data(tree, '.diagnostics.microwave.interferom.phase_comp.ne_l')  # [m^-2]
		(tbeam, vbeam) = get_data(tree, '.oper_diags.ltx_nbi.source_diags.v_hvps')
		tbeamon, tbeamoff = tbeam[np.where(vbeam > 1000)][0], tbeam[np.where(vbeam > 1000)][-1]
		return np.mean(ne[np.where((tt >= tbeamon) & (tt < tbeamoff))])
	except mds.mdsExceptions.TreeNODATA:
		print(f'trouble fetching MDSPlus data for shot {shot}')
		return np.nan


def plot_nbi_rawdata(shots):
	fig, axs = plt.subplots(nrows=4, ncols=2, sharex=True)
	((ax1, ax5), (ax2, ax6), (ax3, ax7), (ax4, ax8)) = axs
	ax1r, ax2r, ax3r, ax4r = ax5.twinx(), ax6.twinx(), ax7.twinx(), ax8.twinx()
	for ax in [ax5, ax6, ax7, ax8]:
		ax.yaxis.set_visible(False)
	ax1.set_ylabel('$V_{arc}$ (V)')
	ax1r.set_ylabel('$I_{arc}$ (A)')
	ax2.set_ylabel('$V_{HVPS}$ (V)')
	ax2r.set_ylabel('$I_{HVPS}$ (A)')
	ax3.set_ylabel('$V_{2Grid}$ (V)')
	ax3r.set_ylabel('$I_{2Grid}$ (A)')
	ax4.set_ylabel('$I_{beam}$ (A)')
	ax4r.set_ylabel('$I_{MIPS}$ (A)')
	ax4.set_xlabel('t (ms)')
	ax8.set_xlabel('t (ms)')
	for axrow in axs:
		for ax in axrow:
			ax.grid()
	for sh in shots:
		try:
			t = get_tree_conn(sh, treename='ltx_nbi')
			print('gathering data for shot {} occurred on {}'.format(sh, get_data(t, '.metadata:timestamp')))
			(tiarc, iarc) = get_data(t, '.source_diags.i_arc')
			(tarc, varc) = get_data(t, '.source_diags.v_arc')
			ax1.plot(tarc * 1.e3, varc, label=sh)
			ax1r.plot(tiarc * 1.e3, iarc)
			(tibeam, ibeam) = get_data(t, '.source_diags.i_hvps')
			(tbeam, vbeam) = get_data(t, '.source_diags.v_hvps')
			ax2.plot(tbeam * 1.e3, vbeam)
			ax2r.plot(tibeam * 1.e3, ibeam)
			(tv2gr, v2gr) = get_data(t, '.source_diags.v_decel_grid')
			(ti2gr, i2gr) = get_data(t, '.source_diags.i_decel_grid')
			ax3.plot(tv2gr * 1.e3, v2gr)
			ax3r.plot(ti2gr * 1.e3, i2gr)
			(timips, imips) = get_data(t, '.source_diags.i_mag_insul')
			(tiacc, iacc) = get_data(t, '.source_diags.i_accel_grid')
			ax4.plot(tiacc * 1.e3, iacc)
			ax4r.plot(timips * 1.e3, imips)
		except:
			print(f'problem encountered processing shot {sh}')
	ax1.legend()
	plt.tight_layout()
	plt.show()


def nbi_ops(shots, nbi_win=None, nbi_tree=False, arc_iv=False, v_thresh=1000.):
	if arc_iv:
		fig, axs = plt.subplots(nrows=4, ncols=2, sharex=True)
		((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = axs
		ax4r = ax4.twinx()
		ax4.yaxis.set_visible(False)
		ax3.set_ylabel('i_arc (A)')
		ax4r.set_ylabel('v_arc (V)')
	else:
		fig, axs = plt.subplots(nrows=3, ncols=2, sharex=True)
		((ax1, ax2), (ax5, ax6), (ax7, ax8)) = axs
	for axrow in axs:
		for ax in axrow:
			ax.grid()
	ax2r = ax2.twinx()
	ax6r = ax6.twinx()
	ax8r = ax8.twinx()
	ax2.yaxis.set_visible(False)
	ax6.yaxis.set_visible(False)
	ax8.yaxis.set_visible(False)
	ax1.set_ylabel('ip (kA)')
	ax2r.set_ylabel('$\int{n_e dl}$ (e18 $m^{-2}$)')
	ax5.set_ylabel('i_beam (A)')
	ax6r.set_ylabel('v_beam (kV)')
	ax7.set_ylabel('perveance (uPerv)')
	ax8r.set_ylabel('power (kW)')
	ax7.set_xlabel('time (s)')
	ax8.set_xlabel('time (s)')
	if nbi_win is not None:
		ax1.set_xlim(nbi_win)
	pav, vav, w = 0, 0, 0
	for sh in shots:
		if sh > 200000:
			tree = 'ltx_nbi'
			nbi_only = True
			prefix = ''
		else:
			tree = 'ltx_b'
			nbi_only = False
			prefix = '.oper_diags.ltx_nbi'
		t = get_tree_conn(sh, treename=tree)
		print('gathering data for shot {} occurred on {}'.format(sh, get_data(t, '.metadata:timestamp')))
		(tibeam, ibeam) = get_data(t, f'{prefix}.source_diags.i_hvps')
		(tbeam, vbeam) = get_data(t, f'{prefix}.source_diags.v_hvps')
		# (tacc, iacc) = get_data(t, f'{prefix}.source_diags.i_accel_grid')
		# (tdec, idec) = get_data(t, f'{prefix}.source_diags.i_decel_grid')
		ibeam = np.interp(tbeam, tibeam, ibeam)
		(tiarc, iarc) = get_data(t, f'{prefix}.source_diags.i_arc')
		(tarc, varc) = get_data(t, f'{prefix}.source_diags.v_arc')
		iarc = np.interp(tarc, tiarc, iarc)
		# iacc = np.interp(tdec, tacc, iacc)
		perv = np.zeros_like(ibeam)
		perv[:] = np.nan
		t_beamon = np.where(vbeam >= v_thresh)  # only look at where beam is above 5kV
		if len(t_beamon[0]) > 0:
			pad = 0.25e-3  # remove this amt from beginning/end of beam
			t_window = np.where((tbeam >= tbeam[t_beamon[0][0]] + pad) & (tbeam <= tbeam[t_beamon[0][-1]] - pad))
			perv[t_window] = ibeam[t_window] / vbeam[t_window] ** 1.5 * 1.e6  # uPerv
		else:
			t_beamon = [.46, .48]  # set some default for setting nbi_win below
		
		if nbi_win is not None:
			twin1, twin2 = nbi_win[0], nbi_win[1]
		else:
			twin1, twin2 = min(t_beamon) - 2.e-3, max(t_beamon) + 2.e-3
		
		# perv = ibeam / vbeam ** 1.5 * 1.e6  # uPerv
		# perv[np.where(vbeam <= v_thresh)] = np.nan
		# else:
		# 	tbeam, ibeam, vbeam, tarc, iarc, varc, perv = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
		pbeam = ibeam * vbeam / 1000.  # kW
		iav = np.where(vbeam > v_thresh)
		vav = (vav * w + np.mean(vbeam[iav])) / (w + 1)
		pav = (pav * w + np.mean(pbeam[iav])) / (w + 1)
		w += 1
		if arc_iv:
			ax3.plot(tarc, iarc)
			ax4r.plot(tarc, varc)
		ax5.plot(tbeam, ibeam)
		# ax5.plot(tdec, iacc, '--')
		# ax5.plot(tdec, iacc-idec, '-.')
		ax6r.plot(tbeam, vbeam / 1000.)
		ax7.plot(tbeam, perv)
		ax8r.plot(tbeam, pbeam)
		if not nbi_only and is_nbi_shot(sh, t):
			(times, ip) = get_data(t, '.diagnostics.magnetic.ip_rog')
			(times2, ne) = get_data(t, '.diagnostics.microwave.interferom.phase_comp.ne_l')  # [m^-2]
			ne = smooth(ne)
			ii = np.where((times >= twin1) & (times <= twin2))
			i2 = np.where((times2 >= twin1) & (times2 <= twin2))
			ax1.plot(times[ii], -ip[ii] / 1000., label=sh)
			ax2r.plot(times2[i2], ne[i2])
	# except:
	# 	print('no tree found for shot {}'.format(sh))
	ax1.legend(ncol=int(np.ceil(len(shots) / 7)))
	# ax2.axis('off')
	
	vav /= 1000.
	print('avg beam voltage: {:.1f} kV'.format(vav))
	ax6r.axhline(vav, c='k', ls='--', alpha=0.5)
	ax6r.annotate('<V>={:.1f}'.format(vav), (twin1, vav))
	print('avg beam power: {:.1f} kW'.format(pav))
	ax8r.axhline(pav, c='k', ls='--', alpha=0.5)
	ax8r.annotate('<P>={:.1f}'.format(pav), (twin1, pav))
	ax1.set_ylim(bottom=0.)
	ax2r.set_ylim(bottom=0.)


# plt.grid(b=True, which='minor')


def plot_nbi(shots):
	nbi_tree_base = '.oper_diags.ltx_nbi.source_diags'
	fig, (ax1, ax2) = plt.subplots(2, 1)
	ax1r = ax1.twinx()
	ax2r = ax2.twinx()
	# ax3r = ax3.twinx()
	ax1.set_ylabel('i_arc')
	ax1r.set_ylabel('v_arc')
	ax2.set_ylabel('i beam')
	ax2r.set_ylabel('v_hvps')
	# ax3.set_ylabel('i_grid')
	# ax3r.set_ylabel('v_decel_grid')
	times = None
	
	def put_data_on_plot(node, ax: Axes, lab=None):
		try:
			(times, dat) = get_data(t, node)
			if lab is None:
				lbl = '{}'.format(sh)
			else:
				lbl = '{}:{}'.format(lab, sh)
			times = times - times[0]
			ax.plot(times, dat, label=lbl)
		except:
			print('no {} data for shot {}'.format(node, sh))
	
	for sh in shots:
		print(sh)
		t = get_tree_conn(sh, treename='ltx_nbi')
		put_data_on_plot(nbi_tree_base + '.i_arc', ax1)
		put_data_on_plot(nbi_tree_base + '.v_arc', ax1r)
		put_data_on_plot(nbi_tree_base + '.i_hvps', ax2)
		put_data_on_plot(nbi_tree_base + '.v_hvps', ax2r)
	# put_data_on_plot(nbi_tree_base + '.i_decel_grid', ax3, lab='decel')
	# put_data_on_plot(nbi_tree_base + '.v_decel_grid', ax3r)
	# put_data_on_plot(nbi_tree_base + '.i_accel_grid', ax3, lab='accel')
	
	ax1.legend()


def arc_v_beam_power_scatter(nbi_shot):
	t = get_tree_conn(nbi_shot)
	nbi_tree_base = '.oper_diags.ltx_nbi.source_diags'
	(t_i_arc, i_arc) = get_data(t, nbi_tree_base + '.i_arc')
	(t_v_arc, v_arc) = get_data(t, nbi_tree_base + '.v_arc')
	(t_i_hvps, i_hvps) = get_data(t, nbi_tree_base + '.i_hvps')
	(t_v_hvps, v_hvps) = get_data(t, nbi_tree_base + '.v_hvps')
	tplt = np.linspace(t_i_arc[0], t_i_arc[-1], num=200)
	i_arc = np.interp(tplt, t_i_arc, i_arc)
	v_arc = np.interp(tplt, t_v_arc, v_arc)
	i_hvps = np.interp(tplt, t_i_hvps, i_hvps)
	v_hvps = np.interp(tplt, t_v_hvps, v_hvps)
	p_arc = i_arc * v_arc * 1.e-6
	p_hvps = i_hvps * v_hvps * 1.e-6
	i_interest = [i for i in range(len(tplt)) if 0.4732 < tplt[i] < 0.4783]
	
	plt.figure('beam power')
	plt.subplot(211)
	plt.plot(tplt, p_arc, label='arc')
	plt.plot(tplt, p_hvps, label='hvps')
	plt.xlabel('time (ms)')
	plt.ylabel('power (MW)')
	plt.legend()
	plt.subplot(212)
	plt.scatter(p_arc, p_hvps)
	plt.scatter(p_arc[i_interest], p_hvps[i_interest], label='beam on')
	plt.xlabel('arc power (MW)')
	plt.ylabel('hvps power (MW)')
	plt.legend()


def ops_scope(shots, nbi_win=False, nbi_tree=False, v_thresh=1000.):
	fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharex=True)
	# ax1r = ax1.twinx()
	# ax2r = ax2.twinx()
	# ax3r = ax3.twinx()
	# ax4r = ax4.twinx()
	ax1.set_ylabel('ip (kA)')
	# ax1r.set_ylabel('v_beam')
	ax3.set_ylabel('i_arc (A)')
	ax4.set_ylabel('v_arc (V)')
	ax5.set_ylabel('i_beam (A)')
	ax6.set_ylabel('v_beam (V)')
	ax7.set_ylabel('perveance (uPerv)')
	ax8.set_ylabel('power (kW)')
	ax7.set_xlabel('time (s)')
	ax8.set_xlabel('time (s)')
	if nbi_win:
		ax1.set_xlim((0.462, 0.469))
	pav, vav, w = 0, 0, 0
	for sh in shots:
		if sh > 200000:
			tree = 'ltx_nbi'
			nbi_only = True
			prefix = ''
		else:
			tree = 'ltx_b'
			nbi_only = False
			prefix = '.oper_diags'
		try:
			t = get_tree_conn(sh, treename=tree)
			print('gathering data for shot {} occurred on {}'.format(sh, get_data(t, '.metadata:timestamp')))
			if not nbi_only and is_nbi_shot(sh, t):
				(times, ip) = get_data(t, '.diagnostics.magnetic.ip_rog')
				ii = np.where((times >= 0.44) & (times <= 0.48))
				ax1.plot(times[ii], -ip[ii] / 1000.)
				ax2.plot(np.nan, np.nan, label=sh)
				(tibeam, ibeam) = get_data(t, '{}.ltx_nbi.source_diags.i_hvps'.format(prefix))
				(tbeam, vbeam) = get_data(t, '{}.ltx_nbi.source_diags.v_hvps'.format(prefix))
				ibeam = np.interp(tbeam, tibeam, ibeam)
				(tiarc, iarc) = get_data(t, '{}.ltx_nbi.source_diags.i_arc'.format(prefix))
				(tarc, varc) = get_data(t, '{}.ltx_nbi.source_diags.v_arc'.format(prefix))
				iarc = np.interp(tarc, tiarc, iarc)
				perv = ibeam / vbeam ** 1.5 * 1.e6  # uPerv
				perv[np.where(vbeam <= v_thresh)] = np.nan
			else:
				tbeam, ibeam, vbeam, tarc, iarc, varc, perv = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
			pbeam = ibeam * vbeam / 1000.  # kW
			iav = np.where(vbeam > v_thresh)
			vav = (vav * w + np.mean(vbeam[iav])) / (w + 1)
			pav = (pav * w + np.mean(pbeam[iav])) / (w + 1)
			w += 1
			ax3.plot(tarc, iarc)
			ax4.plot(tarc, varc)
			ax5.plot(tbeam, ibeam)
			ax6.plot(tbeam, vbeam)
			ax7.plot(tbeam, perv)
			ax8.plot(tbeam, pbeam)
		except:
			print('no tree found for shot {}'.format(sh))
	ax2.legend(ncol=int(np.ceil(len(shots) / 7)))
	ax2.axis('off')
	
	print('avg beam voltage: {} kV'.format(vav / 1000.))
	ax6.axhline(vav, c='k', ls='--', alpha=0.5)
	ax6.annotate('<V>={:.2f}'.format(vav), (np.mean(tbeam), vav))
	print('avg beam power: {} kW'.format(pav))
	ax8.axhline(pav, c='k', ls='--', alpha=0.5)
	ax8.annotate('<P>={:.2f}'.format(pav), (np.mean(tbeam), pav))


def quick_perv(shots):
	for shot in shots:
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
			perv = ones_like(vb)
			perv[:] = np.nan
			perv[t_window] = ib[t_window] / vb[t_window] ** 1.5
			if shot == shots[-1]:
				plt.plot(tv, perv * 1.e6, 'ko-', markevery=25, label=shot)
				plt.axhline(np.nanmean(perv) * 1.e6, ls='--', c='k')
			else:
				plt.plot(tv, perv * 1.e6, label=shot)
		except mds.mdsExceptions.TreeNODATA:
			print(f'trouble fetching MDSPlus data for shot {shot}')
	plt.axhline(15)
	plt.xlabel('time (s)')
	plt.ylabel('perv (e-6)')
	plt.legend()
	plt.ylim(bottom=0)
	plt.show()


def beam_nobeam(shots, twin=[450, 480]):
	fig, ax = plt.subplots()
	axr = ax.twinx()
	for shot in shots:
		try:
			tree = get_tree_conn(shot, treename='ltx_b')
			(tip, ip) = get_data(tree, '.diagnostics.magnetic.ip_rog', times=None)
			(tnel, nel) = get_data(tree, '.diagnostics.microwave.interferom.phase_comp.ne_l')  # [m^-2]
			ax.plot(tip * 1e3, -ip * 1.e-3, label=shot)
			axr.plot(tnel * 1e3, nel)
			if is_nbi_shot(shot, tree):
				(tibeam, ibeam) = get_data(tree, '.oper_diags.ltx_nbi.source_diags.i_hvps')
				(tbeam, vbeam) = get_data(tree, '.oper_diags.ltx_nbi.source_diags.v_hvps')
				ibeam = np.interp(tbeam, tibeam, ibeam)
				ax.plot(tbeam * 1e3, ibeam * vbeam * 1.e-6, 'b')
		except:
			print(f'shot {shot} did not pan out')
	fs = 12
	ax.set_xlim(twin)
	ax.set_ylim((0, 150))
	ax.set_ylabel('$I_p$ (kA)', fontsize=fs)
	ax.set_xlabel('time (ms)', fontsize=fs)
	axr.set_ylabel('$n_eL$ (m$^2$) $P_{nbi}$ (MW)', fontsize=fs)
	axr.set_ylim((0, 1.5e19))
	ax.legend(fontsize=fs)
	ax.tick_params(labelsize=fs)
	axr.tick_params(labelsize=fs)
	plt.tight_layout()
	plt.show()


if __name__ == '__main__':
	"""
	beam on: 461-467
	early = [604] but after TS
	ok = [606,554,555,556,551,553,525,538,542,543,526,527,532,535,536,560,531,544,569,570,571,572,598,599,603]
	bad = [523, 605] but after TS
	high-short = [561,562,563,564,565] CAUTION- using diff beam here! (TS at 463, 464)
	late = [566,567,568,530] CAUTION- beam late! (TS at 465, 466)
	weak-start = [603]
	"""
	# 106525,538,542,543
	nbi_ops(106500 + np.array([25,38,42,43, 36]), nbi_win=[.455, .47], arc_iv=False)
	# nbi_ops(106000 + np.array([606,554,555,556,551,553,525,538,542,543]), nbi_win=[.45, .478], arc_iv=False)
	plt.show()
	
	# for sh in np.arange(509113,509340):
	# 	t = get_tree_conn(sh, treename='ltx_nbi')
	# 	print(f'shot {sh} occurred at {get_data(t, ".metadata:timestamp")}')
	
	# beam_nobeam([103658, 103617])
	# beam_nobeam([103446, 103465])
	# beam_nobeam([103898, 103899])
	# quick_perv(105400 + np.array([28, 29, 30, 31, 32, 33]))
	# nbi_ops([105428, 105427], arc_iv=True, nbi_win=[.45, .475])
	# nbi_ops(105100 + np.array([83, 88, 89]), nbi_win=[.46, .485])
	# nbi_ops(105961+np.arange(6), nbi_win=[.46, .48], arc_iv=True)
	# plot_nbi_rawdata(509000 + np.array([355]))  # 355=fresh Li
	# nbi_ops(106240 + np.array([2, 3, 4, 5, 6, 7, 8]), nbi_win=[.46, .48])
	'''
	nbi_ops([104584], arc_iv=True, nbi_win=[.46,.475])
	nbi_ops([103658, 103617], arc_iv=True)#, nbi_win=[.46, .475])
	# nbi_ops([103465], arc_iv=True, nbi_win=[.46, .475])
	plt.show()
	
	fig, ax = plt.subplots()
	tr = get_tree_conn(101826)
	plot_ip(ax, tr)
	plt.show()
	
	# arc_v_beam_power_scatter(100925)
	get_current_ltx_shot()
	hot_beam = get_shots('hot_notable')
	cold_beam = get_shots('cold_notable')
	hughes_shots = [100981, 100985, 100988]
	# nbi_ops(hughes_shots, nbi_win=(0.466, 0.475))
	nbi_ops(hughes_shots, nbi_win=(0.46, 0.48), arc_iv=True)
'''
