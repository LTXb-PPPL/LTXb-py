import lvm_read
from bills_LTX_MDSplus_toolbox import *
import os

plot_neutralization_fractions = 1
do_11Oct_analysis = 1


def get_calorimeter_dtemp(shot):
	c_ps = 385.  # specific heat of copper J/kg/degC
	m_copp = 1.5834  # calorimeter copper mass [kg]
	treename = 'ltx_b'
	direc = 'Y:/thermocouple/Calorimeter/'
	lvm_files = os.listdir(direc)
	lvm_files = [lvmf for lvmf in lvm_files if lvmf.endswith('.lvm')]
	lvm_times = [datetime.datetime.strptime(lvmf.split('.')[0], '%m%d%Y%H%M%S') for lvmf in lvm_files]
	
	tree = get_tree_conn(shot, treename=treename)
	(ti, ib) = get_data(tree, '.source_diags.i_hvps')
	(tv, vb) = get_data(tree, '.source_diags.v_hvps')
	ts = get_data(tree, '.metadata.timestamp')
	if len(ts.split('/')[0]) < 2:
		ts = f'0{ts}'
	dt = datetime.datetime.strptime(ts, '%m/%d/%Y %I:%M:%S %p')
	sync_file = [lvm_files[i] for i in range(len(lvm_files)) if
	             dt < lvm_times[i] and dt + datetime.timedelta(seconds=60) > lvm_times[i]]
	if len(sync_file) != 1:
		print(f'---- problem syncing shot {shot}, timestamp {ts} with lvm file')
		return np.nan, np.nan
	else:
		print(f'shot {shot} synced with {sync_file[0]}')
		ib = np.interp(tv, ti, ib)  # get ib onto vb axis
		t_beamon = np.where(vb > 5000)  # only look at where beam is above 5kV
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
			return np.nan, np.nan
		else:
			line2fit = (c0[islow_decay] + c1[islow_decay] + c2[islow_decay] + c3[islow_decay] + c4[islow_decay]) / 5.
			time2fit = time[islow_decay]
			fit = np.polyfit(time2fit, line2fit, 1)
			temp_interp = np.interp(time[ipeak], time, time * fit[0] + fit[1])
			dt_meas = temp_interp - baseline
			dt_pred = tot_joules / c_ps / m_copp  # predicted dT (degC) from injected beam energy
			return dt_meas, dt_pred


if plot_neutralization_fractions:
	torr_to_pa = 133.322
	kb = 1.380649e-23  # boltzmann constant [J/K]
	pplot_mtor = np.array([1, 5, 10, 20, 40, 80])  # mTor
	pplot_pa = pplot_mtor * 1.e-3 * torr_to_pa  # [J/m^3]
	tgas = 300.  # [k] room temp neutralizer gas, approx
	
	n = pplot_pa / kb / tgas * 1.e-6  # J/m^3*K/J/K*m^3/cm^3= num/cm^3
	l = 1.  # cm of path in neutralizer
	
	# cross sections from Allison '58
	ebeam = np.array([4, 5, 7, 9, 11, 13, 15, 20, 25, 30, 35])  # keV
	f1_inf = np.array([.095, .095, .095, .095, .1, .105, .125, .185, .220, .280, .320])
	f0_inf = np.array([.895, .891, .886, .885, .88, .875, .855, .797, .764, .706, .680])
	fn1_inf = np.array([.01, .014, .019, .02, .02, .02, .02, .018, .016, .014, 0.])
	sig01 = np.array([4.2, 4.2, 4.2, 4.2, 4.6, 5.0, 5.4, 6.38, 7.28, 8.0, 8.0]) * 1.e-17  # cm^2/atom
	sig10 = np.array([40., 40., 40.2, 40., 39.8, 38.2, 35.8, 29.5, 26., 21., 16.]) * 1.e-17
	ll = np.linspace(1, 100., endpoint=True)
	
	for ie in [4, 7]:  # pick some energies
		fig, ax = plt.subplots()
		for ii, nn in enumerate(n):
			f0nl = f0_inf[ie] * (1. - np.exp(-nn * ll * (sig01[ie] + sig10[ie])))
			ax.plot(ll, f0nl, label=f'{pplot_mtor[ii]}')
		ax.legend(title='neut press [mTor]')
		ax.set_xlabel('neutralizer length (cm)')
		ax.set_ylabel('neutral fraction')
		ax.axvline(55, alpha=0.5, c='k', ls='--')
		ax.set_title(f'Ebeam: {ebeam[ie]}keV')
	a = 1

if do_11Oct_analysis:
	shots = np.arange(104013, 104049)  # shots 104013-104048
	topen = np.array(
		[5, 7.5, 10, 12.5, 15, 2.5, 5, 7.5, 10, 12.5, 15, 2.5, 5, 7.5, 10, 12.5, 15, .1, .2, .4, .8, 1.6, .1, .2, .4,
		 .8, 1.6, 2.5, .2, .4, .4, .4, .4, 1.6, 1.6, 15.])  # puff dt
	toffset = np.zeros_like(topen)
	toffset[0:4] = [10, 7.5, 5, 2.5]  # started puffs at diff times
	psi = np.array(
		[60, 60, 60, 60, 60, 80, 80, 80, 80, 80, 80, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 60, 60, 60,
		 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60])
	i60, i80, i100 = np.where((psi == 60) & (topen > 2))[0], np.where(psi == 80)[0], \
	                 np.where((psi == 100) & (topen > 2))[0]
	i60s, i100s = np.where((psi == 60) & (topen < 2)), (psi == 100) & (topen < 2)
	
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	mig_max = []  # Torr
	for ish, shot in enumerate(shots):
		print(f'{ish}/{len(shots)}')
		tree = get_tree_conn(shot, treename='ltx_b')
		(tig, ig) = get_data(tree, '.oper_diags.ion_gauges.ig_ctrlrs.nbi_tank_mig')
		mig_max.append(max(ig))
	mig_max = np.array(mig_max)
	
	fig, ax1 = plt.subplots()
	f60 = np.polyfit(topen[i60], mig_max[i60] * 1.e3, 1)
	ax1.plot(topen[i60], mig_max[i60] * 1.e3, 'o', label='60psi', c=clrs[0])
	ax1.plot(np.arange(-7, 16), np.arange(-7, 16) * f60[0] + f60[1], '--', c=clrs[0])
	f80 = np.polyfit(topen[i80], mig_max[i80] * 1.e3, 1)
	ax1.plot(topen[i80], mig_max[i80] * 1.e3, 'o', label='80psi', c=clrs[1])
	ax1.plot(np.arange(-7, 16), np.arange(-7, 16) * f80[0] + f80[1], '--', c=clrs[1])
	f100 = np.polyfit(topen[i100], mig_max[i100] * 1.e3, 1)
	ax1.plot(topen[i100], mig_max[i100] * 1.e3, 'o', label='100psi', c=clrs[2])
	ax1.plot(np.arange(-7, 16), np.arange(-7, 16) * f100[0] + f100[1], '--', c=clrs[2])
	ax1.set_xlabel('$t_{open}$ (ms)')
	ax1.set_ylabel('$P_{MIG}$ (mTor)')
	ax1.spines.left.set_position('zero')
	ax1.spines.right.set_color('none')
	ax1.spines.bottom.set_position('zero')
	ax1.spines.top.set_color('none')
	ax1.xaxis.set_ticks_position('bottom')
	ax1.yaxis.set_ticks_position('left')
	plt.legend()
	print(f'60psi y-intercept: {f60[1]:.2f}mTor')
	print(f'80psi intercept: {f80[1]:.2f}mTor')
	print(f'100psi intercept: {f100[1]:.2f}mTor')
	print(f'60psi x-intercept: {-f60[1] / f60[0]:.2f}ms')
	print(f'80psi x-intercept: {-f80[1] / f80[0]:.2f}ms')
	print(f'100psi x-intercept: {-f100[1] / f100[0]:.2f}ms')
	
	plt.figure()
	plt.plot(topen[i60s], mig_max[i60s] * 1.e3, 'o', label='60psi', c=clrs[0])
	plt.plot(np.linspace(0, 1.6, endpoint=True), np.linspace(0, 1.6, endpoint=True) * f60[0] + f60[1], '--', c=clrs[0])
	plt.plot(topen[i100s], mig_max[i100s] * 1.e3, 'o', label='100psi', c=clrs[2])
	plt.plot(np.linspace(0, 1.6, endpoint=True), np.linspace(0, 1.6, endpoint=True) * f100[0] + f100[1], '--',
	         c=clrs[2])
	plt.xlabel('$t_{open}$ (ms)')
	plt.ylabel('$P_{MIG}$ (mTor)')
	plt.legend()
	
	fig, ax = plt.subplots()
	axr = ax.twinx()
	ax.plot([60, 80, 100], [f60[1] * 1.e3 / 60., f80[1] * 1.e3 / 80., f100[1] * 1.e3 / 100.], 'o', c=clrs[0])
	ax.axhline(np.mean([f60[1] * 1.e3 / 60., f80[1] * 1.e3 / 80., f100[1] * 1.e3 / 100.]), c=clrs[0], ls='--')
	axr.plot([60, 80, 100], [f60[0] * 1.e3 / 60., f80[0] * 1.e3 / 80., f100[0] * 1.e3 / 100.], 'o', c=clrs[1])
	axr.axhline(np.mean([f60[0] * 1.e3 / 60., f80[0] * 1.e3 / 80., f100[0] * 1.e3 / 100.]), c=clrs[1], ls='--')
	ax.set_xlabel('bottle pressure (psi)')
	ax.set_ylabel('normalized offset (uTor/psi)', c=clrs[0])
	axr.set_ylabel('normalized slope (mTor/ms/psi)', c=clrs[1])
	ax.set_ylim((0, 4))
	axr.set_ylim((0, .6))
	plt.tight_layout()

plt.show()
