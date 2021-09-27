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


def nbi_ops(shots, nbi_win=None, nbi_tree=False, arc_iv=False, v_thresh=1000.):
	if arc_iv:
		fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharex=True)
		ax4r = ax4.twinx()
		ax4.yaxis.set_visible(False)
		ax3.set_ylabel('i_arc (A)')
		ax4r.set_ylabel('v_arc (V)')
	else:
		fig, ((ax1, ax2), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=3, ncols=2, sharex=True)
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
			prefix = '.oper_diags'
		try:
			t = get_tree_conn(sh, treename=tree)
			print('gathering data for shot {} occurred on {}'.format(sh, get_data(t, '.metadata:timestamp')))
			if not nbi_only and is_nbi_shot(sh, t):
				(times, ip) = get_data(t, '.diagnostics.magnetic.ip_rog')
				(times2, ne) = get_data(t, '.diagnostics.microwave.interferom.phase_comp.ne_l')  # [m^-2]
				ne = smooth(ne)
				if nbi_win is not None:
					twin1, twin2 = nbi_win[0], nbi_win[1]
				else:
					twin1, twin2 = min(times), max(times)
				ii = np.where((times >= twin1) & (times <= twin2))
				i2 = np.where((times2 >= twin1) & (times2 <= twin2))
				ax1.plot(times[ii], -ip[ii] / 1000., label=sh)
				ax2r.plot(times2[i2], ne[i2])
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
			if arc_iv:
				ax3.plot(tarc, iarc)
				ax4r.plot(tarc, varc)
			ax5.plot(tbeam, ibeam)
			ax6r.plot(tbeam, vbeam / 1000.)
			ax7.plot(tbeam, perv)
			ax8r.plot(tbeam, pbeam)
		except:
			print('no tree found for shot {}'.format(sh))
	ax1.legend(ncol=int(np.ceil(len(shots) / 7)))
	# ax2.axis('off')
	
	vav /= 1000.
	print('avg beam voltage: {} kV'.format(vav))
	ax6r.axhline(vav, c='k', ls='--', alpha=0.5)
	ax6r.annotate('<V>={:.2f}'.format(vav), (twin1, vav))
	print('avg beam power: {} kW'.format(pav))
	ax8r.axhline(pav, c='k', ls='--', alpha=0.5)
	ax8r.annotate('<P>={:.2f}'.format(pav), (twin1, pav))
	ax1.set_ylim(bottom=0.)
	ax2r.set_ylim(bottom=0.)


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


if __name__ == '__main__':
	nbi_ops([102796], arc_iv=True)
	nbi_ops([103617], arc_iv=True)#, nbi_win=[.46, .475])
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
	plt.show()
