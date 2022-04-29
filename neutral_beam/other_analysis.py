import glob

import matplotlib.pyplot as plt
import numpy as np
from helpful_stuff import is_nbi_shot, is_good_shot
from bills_LTX_MDSplus_toolbox import get_tree_conn, get_data
import os
import plotly
import json
from helpful_stuff import read_eqdsk2
from conbeam import conbeam_tools
from scipy.io import loadmat


def rmag_v_rlcfs(rescan=True, ploteach=False):
	# look through equilibria, plot r_mag vs r_lcfs with eq identification (plotly?)
	# then we can analyze eq with same r_mag and diff r_lcfs and vice versa to interrogate coupling dependencies
	# look at stuff in ltx_beambot if you only want eqs during NBI for some reason
	
	# rescanning takes a while, so we'll do it if necessary but save the data for quick reruns
	fsav = 'Z:/PycharmProjects/LTXb-py/neutral_beam/rmag_v_rlcfs.npz'
	if rescan:
		reco_dir = 'Y:/reconstructions/runday/'
		
		print(f'gathering eqdsks...')
		eqdsks = glob.glob(f'{reco_dir}**/LTX*.eqdsk', recursive=True)
		print(f'...complete')
		# eqdsks = glob.glob(f'{reco_dir}SHOTS/105894/**/LTX*.eqdsk', recursive=True)  # fast scan for debugging
		rmag, rlcfs = [], []
		for i, eq in enumerate(eqdsks):
			eq = eq.replace('\\', '/')
			print(f'extracting info from {eq}, #{i}/{len(eqdsks)}')
			out = read_eqdsk2(eq)
			if ploteach:
				r_1d = np.linspace(out['rleft'], out['rleft'] + out['rdim'], num=out['nr'])
				z_1d = np.linspace(out['zmid'] - out['zdim'] / 2, out['zmid'] + out['zdim'] / 2, num=out['nz'])
				R, Z = np.meshgrid(r_1d, z_1d)
				plt.contour(R, Z, out['psirz'].transpose())
				plt.plot(out['raxis'], out['zaxis'], 'ks')
				a = 1
			rmag.append(out['raxis'])
			rlcfs.append(np.max(out['rzout'][0, :]))
		print(f'--saved data--')
		np.savez(fsav, rmag=rmag, rlcfs=rlcfs, eqdsks=eqdsks)
		rmag, rlcfs, eqdsks = np.array(rmag), np.array(rlcfs), np.array(eqdsks)
	else:
		dat = np.load(fsav)
		rmag = np.array(dat['rmag'])
		rlcfs = np.array(dat['rlcfs'])
		eqdsks = np.array(dat['eqdsks'])
	
	tol = .001  # 1 cm
	rmagchoice = .4
	rlcfschoice = .57
	irmag = np.where(abs(rmag - rmagchoice) < tol)[0]
	irlcfs = np.where(abs(rlcfs - rlcfschoice) < tol)[0]
	print(f'found {len(irmag)} pts within {tol} of rmag')
	print(f'found {len(irlcfs)} pts within {tol} of rlcfs')
	
	fig, ax = plt.subplots()
	ax.plot(rlcfs, rmag, 'o', zorder=1)
	ax.fill_between(ax.get_xlim(), [rmagchoice - tol, rmagchoice - tol], [rmagchoice + tol, rmagchoice + tol],
	                facecolor='k', alpha=.25, zorder=2)
	yl = ax.get_ylim()
	ax.fill_between([rlcfschoice - tol, rlcfschoice + tol], [yl[0], yl[0]], [yl[1], yl[1]], facecolor='k', alpha=.25,
	                zorder=2)
	ax.set_ylim(yl)
	ax.set_xlabel('$r_{lcfs}$')
	ax.set_ylabel('$r_{mag}$')
	
	fig2, (ax1, ax2) = plt.subplots(ncols=2, sharey=True)
	eqs = eqdsks[irmag]
	rlcfs_pts = rlcfs[irmag]
	runnums = [int(eq.split('_')[-1].split('.')[0]) for eq in eqs]
	shots = [int(eq.split('_')[-2]) for eq in eqs]
	shine, prompt, coup = [], [], []
	for sh, rn, eq in zip(shots, runnums, eqs):
		matfile = f'Z:/matlab/botdata/{sh}_{rn:02}.mat'
		if not os.path.exists(matfile):
			conbeam_tools.create_map(eq, sh, rn)
		try:
			matdat = loadmat(matfile)
			shine.append(matdat['shinethru'][0][0])
			prompt.append(matdat['prompt_loss'][0][0])
			coup.append(matdat['coupled'][0][0])
		except FileNotFoundError:
			shine.append(np.nan)
			prompt.append(np.nan)
			coup.append(np.nan)
	ax1.plot(rlcfs_pts, shine, label='shinethrough')
	ax1.plot(rlcfs_pts, prompt, label='prompt loss')
	ax1.plot(rlcfs_pts, coup, label='coupled')
	ax1.set_title(f'rmag = {rmagchoice}+/-{tol}')

	eqs = eqdsks[irlcfs]
	rmag_pts = rmag[irlcfs]
	runnums = [int(eq.split('_')[-1].split('.')[0]) for eq in eqs]
	shots = [int(eq.split('_')[-2]) for eq in eqs]
	shine, prompt, coup = [], [], []
	for sh, rn, eq in zip(shots, runnums, eqs):
		matfile = f'Z:/matlab/botdata/{sh}_{rn:02}.mat'
		if not os.path.exists(matfile):
			conbeam_tools.create_map(eq, sh, rn)
		try:
			matdat = loadmat(matfile)
			shine.append(matdat['shinethru'][0][0])
			prompt.append(matdat['prompt_loss'][0][0])
			coup.append(matdat['coupled'][0][0])
		except FileNotFoundError:
			shine.append(np.nan)
			prompt.append(np.nan)
			coup.append(np.nan)
	ax2.plot(rmag_pts, shine, label='shinethrough')
	ax2.plot(rmag_pts, prompt, label='prompt loss')
	ax2.plot(rmag_pts, coup, label='coupled')
	ax2.set_title(f'rlcfs = {rlcfschoice}+/-{tol}')

	ax1.legend()
	ax1.set_xlabel('$R_{lcfs}$ (m)')
	ax1.set_ylabel('beam fraction')
	ax2.set_xlabel('$R_{mag}$ (m)')
	plt.tight_layout()
	
	plt.show()


def density_analysis_8feb22():
	shots_8feb = np.arange(25) + 105177  # shots 105177-105201
	find_good_nbi = 0
	if find_good_nbi:
		good_nbi = []
		for shot in shots_8feb:
			t = get_tree_conn(shot, treename='ltx_b')
			if is_nbi_shot(shot, t):
				good_nbi.append(shot)
	else:
		good_nbi = [105177, 105178, 105179, 105180, 105181, 105182, 105183, 105184, 105185, 105186, 105187, 105188,
		            105189, 105190, 105191, 105192, 105193, 105194, 105195, 105196, 105197, 105198, 105199, 105200,
		            105201]
	fig, (ax1, ax2) = plt.subplots(ncols=2)
	ax2r = ax2.twinx()
	for shot in good_nbi:
		print(shot)
		try:
			t = get_tree_conn(shot, treename='ltx_b')
			# (tip, ip) = get_data(t, '.diagnostics.magnetic.ip_rog')
			(tnel, nel) = get_data(t, '.diagnostics.microwave.interferom.phase_comp.ne_l')  # [m^-2]
			(tibeam, ibeam) = get_data(t, '.oper_diags.ltx_nbi.source_diags.i_hvps')
			(tbeam, vbeam) = get_data(t, '.oper_diags.ltx_nbi.source_diags.v_hvps')
			ibeam = np.interp(tbeam, tibeam, ibeam)
			pbeam = ibeam * vbeam
			ax1.plot(tbeam, pbeam)
			ax2.plot(tnel, nel)
		# ax2r.plot(tip, ip)
		except:
			print(f'problem with {shot}')
	for ax in [ax1, ax2]:
		ax.set_xlim((.44, .5))
	plt.tight_layout()
	plt.show()


if __name__ == '__main__':
	# density_analysis_8feb22()
	print('running main')
	rmag_v_rlcfs(rescan=False)
