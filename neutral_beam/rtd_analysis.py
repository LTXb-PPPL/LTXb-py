import matplotlib.pyplot as plt
import numpy as np
from helpful_stuff import SimpleSignal
import lvm_read
from bills_LTX_MDSplus_toolbox import *
import pickle

# shots = [103751, 103752, 103753]
shots_since_9_1_21 = [103734, 103735, 103736, 103737, 103738, 103739, 103740, 103742, 103743, 103744, 103745, 103746,
                      103747, 103748, 103749, 103750, 103751, 103752, 103753, 103754, 103756, 103757, 103758, 103759,
                      103775, 103776, 103777, 103778, 103779, 103780, 103781, 103782, 103783, 103784, 103794, 103797,
                      103798, 103799, 103800, 103801, 103802, 103838, 103839, 103840, 103841, 103842, 103843, 103844,
                      103845, 103847, 103848, 103849, 103850, 103851, 103853, 103854, 103855, 103856, 103857, 103859,
                      103860, 103873, 103874, 103875, 103876, 103877, 103878, 103879, 103884, 103885, 103886, 103887,
                      103888, 103889, 103890, 103891, 103892, 103894, 103895, 103896, 103897, 103898, 103905, 103906,
                      103907, 103911, 103912, 103913, 103914, 103915, 103916, 103917, 103918, 103919, 103920, 103921,
                      103922]

shots_since_8_1_21 = [103879, 103878, 103877, 103876, 103875, 103874, 103873, 103860, 103859, 103857, 103856, 103855,
                      103854, 103853, 103851, 103850, 103849, 103848, 103847, 103845, 103844, 103843, 103842, 103841,
                      103840, 103839, 103838, 103802, 103801, 103800, 103799, 103798, 103797, 103794, 103784, 103783,
                      103782, 103781, 103780, 103779, 103778, 103777, 103776, 103775, 103759, 103758, 103757, 103756,
                      103754, 103753, 103752, 103751, 103750, 103749, 103748, 103747, 103746, 103745, 103744, 103743,
                      103742, 103740, 103739, 103738, 103737, 103736, 103735, 103734, 103637, 103636, 103635, 103634,
                      103633, 103632, 103631, 103630, 103629, 103628, 103627, 103626, 103625, 103624, 103623, 103622,
                      103621, 103620, 103619, 103618, 103617, 103616, 103615, 103614, 103613, 103612, 103611, 103610,
                      103609, 103608, 103607, 103606, 103605, 103604, 103603, 103602, 103595, 103593, 103592, 103591,
                      103590, 103589, 103588, 103587, 103586, 103585, 103583, 103582, 103581, 103580, 103579, 103578,
                      103577, 103576, 103575, 103572, 103571, 103570, 103569, 103563, 103562, 103561, 103560]
shots = shots_since_8_1_21
dir = 'Y:/RTD/'
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']

dbug = 0
if dbug:
	fig2, dax = plt.subplots()


def do_rtd_analysis():
	oldshots = 104400 + np.array([18, 20, 21])
	newshots = 104400 + np.array([15, 16, 17])
	
	def rtd_analysis(shots, m, legend=False):
		norm_rtd = np.zeros((len(shots), 19))  # index of [shot, rtd]
		rtd_arr = np.zeros((len(shots), 19))
		ja, pa = [], []  # shot arrays for joules, av perv
		global ax1, ax2, ax3, fig2, ax11, ax22, ax33
		for ish, shot in enumerate(shots):
			lvm_file = f'{dir}{shot}.lvm'
			if os.path.exists(lvm_file):
				print(shot)
				tree = get_tree_conn(shot, treename='ltx_b')
				(ti, ib) = get_data(tree, '.oper_diags.ltx_nbi.source_diags.i_hvps')
				(tv, vb) = get_data(tree, '.oper_diags.ltx_nbi.source_diags.v_hvps')
				ib = np.interp(tv, ti, ib)  # get ib onto vb axis
				t_beamon = np.where(vb > 5000)  # only look at where beam is above 5kV
				if len(t_beamon[0]) < 10:
					ja.append(np.nan)
					pa.append(np.nan)
					norm_rtd[ish, :] = np.nan
					rtd_arr[ish, :] = np.nan
				else:
					perv = ones_like(vb)
					pb = ib * vb  # beam power [W]
					perv[:] = np.nan
					perv[t_beamon] = ib[t_beamon] / vb[t_beamon] ** 1.5
					av_perv = np.mean(perv[np.where((tv >= tv[t_beamon[0][0]]) & (tv <= tv[t_beamon[0][-1]]))])
					tot_joules = np.sum([pb[i] * (tv[i] - tv[i - 1]) for i in np.arange(1, len(tv))])
					ja.append(tot_joules)
					pa.append(av_perv * 1.e6)
					lvm = lvm_read.read(lvm_file, dump_file=False, read_from_pickle=False)
					lvm0 = lvm[0]
					rtd_nam = lvm0['Channel names'][1:20]  # u1-4,l1-4,d1-11
					for i in np.arange(1, 20):
						norm_rtd[ish, i - 1] = (max(lvm0['data'][:, i]) - min(lvm0['data'][:, i])) / tot_joules
						rtd_arr[ish, i - 1] = (max(lvm0['data'][:, i]) - min(lvm0['data'][:, i]))
						if dbug:
							dax.plot(lvm0['data'][:, 8:11], label=rtd_nam[8:11])
			else:
				ja.append(np.nan)
				pa.append(np.nan)
				norm_rtd[ish, :] = np.nan
				rtd_arr[ish, :] = np.nan
		# bad signals on D1, D3
		# norm_rtd[:, 8] = np.nan
		# norm_rtd[:, 10] = np.nan
		# rtd_arr[:, 8] = np.nan
		# rtd_arr[:, 10] = np.nan
		
		# for ax in [ax1, ax2, ax3]:
		# 	ax.set_xlim(left=7.5)
		# ax.set_ylim(top=2.)
		if legend:
			ax1.plot(pa, norm_rtd[:, 0:4] * 1.e3, m, label=rtd_nam[0:4])  # upper scrapers
			ax2.plot(pa, norm_rtd[:, 4:8] * 1.e3, m, label=rtd_nam[4:8])  # lower scrapers
			ax3.plot(pa, norm_rtd[:, 8:19] * 1.e3, m, label=rtd_nam[8:19])  # dump
			plt.tight_layout()
			ax11.plot(ja, rtd_arr[:, 0:4], m, label=rtd_nam[0:4])  # upper scrapers
			ax22.plot(ja, rtd_arr[:, 4:8], m, label=rtd_nam[4:8])  # lower scrapers
			ax33.plot(ja, rtd_arr[:, 8:19], m, label=rtd_nam[8:19])  # dump
		else:
			ax1.plot(pa, norm_rtd[:, 0:4] * 1.e3, m)  # upper scrapers
			ax2.plot(pa, norm_rtd[:, 4:8] * 1.e3, m)  # lower scrapers
			ax3.plot(pa, norm_rtd[:, 8:19] * 1.e3, m)  # dump
			plt.tight_layout()
			ax11.plot(ja, rtd_arr[:, 0:4], m)  # upper scrapers
			ax22.plot(ja, rtd_arr[:, 4:8], m)  # lower scrapers
			ax33.plot(ja, rtd_arr[:, 8:19], m)  # dump
		return ja, rtd_arr, pa, norm_rtd, rtd_nam
	
	fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(12, 4))
	fig2, (ax11, ax22, ax33) = plt.subplots(ncols=3, sharey=True, figsize=(12, 4))
	oja, ortd_arr, opa, onorm_rtd, rtd_nam = rtd_analysis(oldshots, 'o', legend=True)
	for ax in [ax1, ax2, ax3, ax11, ax22, ax33]:
		ax.set_prop_cycle(None)  # reset colors
	nja, nrtd_arr, npa, nnorm_rtd, _ = rtd_analysis(newshots, 's')
	
	# ax1t = ax1.twiny()
	# ax1t.plot(pa, np.transpose(np.transpose(np.outer(np.ones(3), rtd_arr[2, 0:4])) * ja / ja[2]), '--')
	# ax1t.set_xlabel('Perv (e-6)')
	# ax2t = ax2.twiny()
	# ax2t.plot(pa, np.transpose(np.transpose(np.outer(np.ones(3), rtd_arr[2, 4:8])) * ja / ja[2]), '--')
	# ax2t.set_xlabel('Perv (e-6)')
	# ax3t = ax3.twiny()
	# ax3t.plot(pa, np.transpose(np.transpose(np.outer(np.ones(3), rtd_arr[2, 8:19])) * ja / ja[2]), '--')
	# ax3t.set_xlabel('Perv (e-6)')
	
	leg_fs = 10
	fs = leg_fs + 2
	ax1.legend(ncol=2, fontsize=leg_fs)
	ax2.legend(ncol=2, fontsize=leg_fs)
	ax3.legend(ncol=3, fontsize=leg_fs)
	ax11.legend(ncol=2, fontsize=leg_fs)
	ax22.legend(ncol=2, fontsize=leg_fs)
	ax33.legend(ncol=3, fontsize=leg_fs)
	for ax in [ax1, ax2, ax3, ax11, ax22, ax33]:
		ax.tick_params(labelsize=fs)
	
	ax1.set_ylabel('$\Delta T/E_{tot}$ ($\degree$C/kJ)', fontsize=fs)
	ax2.set_xlabel('Perv (e-6)', fontsize=fs)
	ax1.set_title('Upper', fontsize=fs)
	ax2.set_title('Lower', fontsize=fs)
	ax3.set_title('Dump', fontsize=fs)
	ax11.set_ylabel('$\Delta T$ ($\degree$C)', fontsize=fs)
	ax22.set_xlabel('$E_{tot}$ (J)', fontsize=fs)
	ax11.set_title('Upper', fontsize=fs)
	ax22.set_title('Lower', fontsize=fs)
	ax33.set_title('Dump', fontsize=fs)
	
	plt.tight_layout()
	
	ortda = np.mean(ortd_arr, axis=0)
	nrtda = np.mean(nrtd_arr, axis=0)
	fig3, axx = plt.subplots()
	fs = fs + 2
	axx.plot(rtd_nam, nrtda / ortda, 'o')
	axx.set_ylabel('Temp Increase Factor', fontsize=fs)
	axx.tick_params(labelsize=fs)
	axx.set_xticklabels(rtd_nam, rotation=90)


def plot_rtd_sigs(shot):
	lvm_file = f'{dir}{shot}.lvm'
	if os.path.exists(lvm_file):
		print(shot)
		# tree = get_tree_conn(shot, treename='ltx_b')
		# (ti, ib) = get_data(tree, '.oper_diags.ltx_nbi.source_diags.i_hvps')
		# (tv, vb) = get_data(tree, '.oper_diags.ltx_nbi.source_diags.v_hvps')
		# ib = np.interp(tv, ti, ib)  # get ib onto vb axis
		# t_beamon = np.where(vb > 5000)  # only look at where beam is above 5kV
		# perv = ones_like(vb)
		# pb = ib * vb  # beam power [W]
		# perv[:] = np.nan
		# perv[t_beamon] = ib[t_beamon] / vb[t_beamon] ** 1.5
		# av_perv = np.mean(perv[np.where((tv >= tv[t_beamon[0][0]]) & (tv <= tv[t_beamon[0][-1]]))])
		# tot_joules = np.sum([pb[i] * (tv[i] - tv[i - 1]) for i in np.arange(1, len(tv))])
		lvm = lvm_read.read(lvm_file, dump_file=False, read_from_pickle=False)
		lvm0 = lvm[0]
		rtd_nam = lvm0['Channel names'][1:20]  # u1-4,l1-4,d1-11
		for i in np.arange(1, 9):
			plt.plot(lvm0['data'][:, 0], lvm0['data'][:, i], label=rtd_nam[i - 1])
		plt.legend()
		plt.xlabel('Time (s)')
		plt.ylabel('Temp (degC)')


def avg_perv(shot, Ej=False):
	tree = get_tree_conn(shot, treename='ltx_b')
	(ti, ib) = get_data(tree, '.oper_diags.ltx_nbi.source_diags.i_hvps')
	(tv, vb) = get_data(tree, '.oper_diags.ltx_nbi.source_diags.v_hvps')
	ib = np.interp(tv, ti, ib)  # get ib onto vb axis
	# t_beamon = np.where(vb > 5000)  # only look at where beam is above 5kV
	t_beamon = np.where((.46 < tv) & (tv < .463))  # look between 460-463ms for these shots (ignore rough Ibeam startup)
	pb = ib * vb  # beam power [W]
	tot_joules = np.sum([pb[i] * (tv[i] - tv[i - 1]) for i in np.arange(1, len(tv))])
	if len(t_beamon[0]) < 10:
		av_perv = np.nan
	else:
		perv = ones_like(vb)
		perv[:] = np.nan
		perv[t_beamon] = ib[t_beamon] / vb[t_beamon] ** 1.5
		av_perv = np.nanmean(perv)
	if Ej:
		return av_perv, tot_joules
	else:
		return av_perv


def get_rtd_quick_response(shots):
	tp = np.linspace(-20, 100, num=1000, endpoint=True)
	num = 0
	up = np.zeros((len(tp), 4))
	dn = np.zeros((len(tp), 4))
	for ish, shot in enumerate(shots):
		_, tot_j = avg_perv(shot, Ej=True)
		lvm_file = f'{dir}{shot}.lvm'
		if os.path.exists(lvm_file):
			print(shot)
			lvm = lvm_read.read(lvm_file, dump_file=False, read_from_pickle=False)
			lvm0 = lvm[0]
			rtd_nam = lvm0['Channel names'][1:20]  # u1-4,l1-4,d1-11
			# sigs 2,3,4 on uppers seem good to use to locate tbeam0
			ibeam = np.where(np.diff(lvm0['data'][:, 2]) == max(np.diff(lvm0['data'][:, 2])))[0][0]
			t = lvm0['data'][:, 0] - lvm0['data'][:, 0][ibeam]
			for iscrp in [0, 1, 2, 3]:
				up[:, iscrp] += np.interp(tp, t, lvm0['data'][:, iscrp + 1] - mean(
					lvm0['data'][:ibeam - 5, iscrp + 1])) / tot_j
				dn[:, iscrp] += np.interp(tp, t, lvm0['data'][:, iscrp + 5] - mean(
					lvm0['data'][:ibeam - 5, iscrp + 5])) / tot_j
			num += 1
	up /= num
	dn /= num
	
	# sigs 2,3,4 on uppers seem good to use to locate tbeam0
	ibeam = np.where(np.diff(up[:, 1]) == max(np.diff(up[:, 1])))[0][0]  # use scraper number 2
	quick = 10  # num of seconds to wait for temp response
	iquick = int(quick / (tp[1] - tp[0]))  # num of data points in quick
	dt = np.zeros(8)
	for i in range(4):
		dt[i] = max(up[ibeam:ibeam + iquick, i])
		dt[i + 4] = max(dn[ibeam:ibeam + iquick, i])
	return tp, up, dn, None, dt


def perveance_scan_rtd_analysis(shots):
	# plot temp rise on beam dump/scrapers vs avg perveance
	rtd_arr = np.zeros((19, len(shots)))
	av_pervs = np.zeros(len(shots))
	for ish, shot in enumerate(shots):
		lvm_file = f'{dir}{shot}.lvm'
		if os.path.exists(lvm_file):
			print(shot)
			avp, totj = avg_perv(shot, Ej=True)
			av_pervs[ish] = avp
			lvm = lvm_read.read(lvm_file, dump_file=False, read_from_pickle=False)
			lvm0 = lvm[0]
			rtd_nam = lvm0['Channel names'][1:20]  # u1-4,l1-4,d1-11
			for i in np.arange(1, 20):
				rtd_arr[i - 1, ish] = (max(lvm0['data'][:, i]) - min(lvm0['data'][:, i])) / totj
		else:
			av_pervs[ish] = np.nan
			rtd_arr[:, ish] = np.nan
	
	isort = np.argsort(av_pervs)
	fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True)
	for irtd in [0, 1, 2, 3]:  #
		ax1.plot(av_pervs[isort] * 1.e6, rtd_arr[irtd, :][isort] * 1.e3, 'o-', label=rtd_nam[irtd])
	for irtd in [4, 5, 6, 7]:
		ax2.plot(av_pervs[isort] * 1.e6, rtd_arr[irtd, :][isort] * 1.e3, 'o-', label=rtd_nam[irtd])
	for irtd in [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]:
		ax3.plot(av_pervs[isort] * 1.e6, rtd_arr[irtd, :][isort] * 1.e3, 'o-', label=rtd_nam[irtd])
	for ax in [ax1, ax2, ax3]:
		ax.legend()
	ax1.set_title('Upper Scrapers')
	ax2.set_title('Lower Scrapers')
	ax3.set_title('Far Side Dump')
	ax2.set_xlabel('Perveance ($x10^{-6}$)')
	ax1.set_ylabel('$\Delta T/E_{tot}$ ($\degree$C/kJ)')
	
	fig2, ax4 = plt.subplots()
	ax4.plot(av_pervs[isort] * 1.e6, rtd_arr[0, :][isort] / rtd_arr[4, :][isort], 'o-', label='1')
	ax4.plot(av_pervs[isort] * 1.e6, rtd_arr[1, :][isort] / rtd_arr[5, :][isort], 'o-', label='2')
	ax4.plot(av_pervs[isort] * 1.e6, rtd_arr[2, :][isort] / rtd_arr[6, :][isort], 'o-', label='3')
	ax4.plot(av_pervs[isort] * 1.e6, rtd_arr[3, :][isort] / rtd_arr[7, :][isort], 'o-', label='4')
	ax4.set_xlabel('Perveance ($x10^{-6}$)')
	ax4.set_ylabel('Upper/Lower ratio')
	ax4.set_ylim(bottom=0)
	ax4.legend(title='RTD #')


def compare_rtd_signals():
	import matplotlib.gridspec as gridspec
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	
	ir_24jan22 = 104900 + np.array([6, 8, 9, 11, 13, 17, 18, 20, 23])
	ir_26jan22 = np.arange(15) + 104960
	
	tp, preup, predn, predmp, dt_pre = get_rtd_quick_response(ir_24jan22)
	tp, postup, postdn, postdmp, dt_post = get_rtd_quick_response(ir_26jan22)
	
	fig = plt.figure(tight_layout=True, figsize=(10, 10))
	gs = gridspec.GridSpec(3, 2)
	ax1 = fig.add_subplot(gs[0, 0])
	ax2 = fig.add_subplot(gs[0, 1])
	ax3 = fig.add_subplot(gs[1, 0])
	ax4 = fig.add_subplot(gs[1, 1])
	ax5 = fig.add_subplot(gs[2, :])
	for i in [0, 1, 2, 3]:
		ax1.plot(tp, preup[:, i] * 1.e3, clrs[i], label=str(i))
		ax2.plot(tp, postup[:, i] * 1.e3, clrs[i])
		ax3.plot(tp, predn[:, i] * 1.e3, clrs[i+4], label=str(i))
		ax4.plot(tp, postdn[:, i] * 1.e3, clrs[i+4])
	ax1.legend(title='scraper #')
	ax3.legend(title='scraper #')
	
	ax5.plot([4, 3, 2, 1], dt_post[:4] * 1.e3 - dt_pre[:4] * 1.e3, 's-', label='upper')
	ax5.plot([4, 3, 2, 1], dt_post[4:] * 1.e3 - dt_pre[4:] * 1.e3, 's-', label='lower')
	ax5.legend()
	
	ax1.set_title('Before')
	ax2.set_title('After')
	ax1.set_ylabel('UPPER (deg C per kJ)')
	ax1.set_xlabel('time rel to beam')
	ax1.set_ylim(ax2.get_ylim())
	ax2.set_xlabel('time rel to beam')
	ax3.set_ylabel('LOWER (deg C per kJ)')
	ax3.set_xlabel('time rel to beam')
	ax4.set_ylim(ax3.get_ylim())
	ax4.set_xlabel('time rel to beam')
	ax5.set_ylabel('Change in dT due to realignment\n(deg C per kJ)')
	ax5.set_xticks([1, 2, 3, 4])
	ax5.set_xticklabels(['4', '3', '2', '1'])
	ax5.set_xlabel('Scraper RTD number')
	plt.show()


if __name__ == '__main__':
	# do_rtd_analysis()
	# plot_rtd_sigs(104849)
	perv_scan_14Jan22 = 104800 + np.array([18, 22, 26, 30, 31, 32, 33, 37, 38, 39, 42, 44, 49, 50, 53, 54, 56, 57, 58,
	                                       59, 61, 62, 63, 64, 65])
	# perveance_scan_rtd_analysis(perv_scan_14Jan22)
	
	ir_24jan22 = 104900 + np.array([6, 8, 9, 11, 13, 17, 18, 20, 23])
	# perveance_scan_rtd_analysis(ir_24jan22)
	ir_26jan22 = np.arange(15) + 104960
	# perveance_scan_rtd_analysis(ir_26jan22)
	
	compare_rtd_signals()
	plt.show()
