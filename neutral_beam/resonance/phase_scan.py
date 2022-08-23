import matplotlib.pyplot as plt
import numpy as np
import glob
import os

from helpful_stuff import closest
from neutral_beam.resonance.resonance_tools import convert_sav_to_dict
import pickle
from neutral_beam.resonance.vdote import vdote

if os.path.exists('Z:/users/wcapecch/'):
	direc = 'Z:/users/wcapecch/'
	proj_direc = 'Z:/PycharmProjects/LTXb-py/'
	pych_direc = 'Z:/PycharmProjects/'
elif os.path.exists('//samba/wcapecch/'):
	direc = '//samba/wcapecch/'
	proj_direc = '//samba/wcapecch/PycharmProjects/LTXb-py/'
	pych_direc = '//samba/wcapecch/PycharmProjects/'


def phase_scan_ltx(mnum=1, nnum=2, fepm=None, psimode=None, phi0=None, verbose=False):  # , parallel=parallel,
	# ps=ps, click=click, newdata=newdata, contours=contours, thesis=thesis, reclick=reclick)
	
	nf = 51
	nr = 50
	no = 1
	
	if psimode is None:
		psimode = np.linspace(0, 1., num=nr)
	if fepm is None:
		fepm = np.linspace(1.e3, 100.e3, num=nf)
	if phi0 is None:
		phi0 = [0.]
	
	orbit_dir = f'{pych_direc}data/orbits/'
	orbits = glob.glob(f'{orbit_dir}*orbit.sav')
	for orbit in orbits:
		if not os.path.isfile(f'{orbit[:-4]}.pkl'):  # pickled version doesn't exist
			convert_sav_to_dict(orbit)
	orbits = glob.glob(f'{orbit_dir}*orbit.pkl')  # list of pickled orbits
	
	psi_res = []  # array of computed resonance psi radii
	f_res = []  # same for freqs of epm
	o_res = []  # same for phi0 of epm
	# r_orb = []  # array of orbit computed radii
	f_orb = np.zeros((3, len(orbits)))  # orbit predicted fepm for l=-1,0,1
	ioniz0 = np.zeros((4, len(orbits)))  # x,y,R,z of ionization pt for each orbit
	# gc0 = np.zeros((4, len(orbits)))  # x,y,R,z of guiding center pt for each orbit at ionization
	
	for o in range(len(orbits)):
		orbit_dat = pickle.load(open(orbits[o], 'rb'))
		orb = orbit_dat['orb']
		ioniz0[:, o] = [orb['X'][0], orb['Y'][0], orb['R'][0], orb['Z'][0]]
		# TODO compute GC x,y,R,z and update gc0 array- DO I NEED GC DATA? maybe for plotting initial GC at ionization?
		print(orbits[o])
		if 'resonance' not in orbit_dat.keys():
			print(f'{len(orbits) - o} orbits remaining')
			recon = orbit_dat['orb']['RECONST_FILE'].decode()
			if not os.path.exists(recon):
				eqdsk_shrt = recon.split('/')[-1].split('\\')[-1]
				recon = f'{direc}/datasets/{eqdsk_shrt}'
				if not os.path.exists(recon):
					print(f"can't find reconstruction file: {recon}")
					return
			resonance = vdote(mnum, nnum, fepm, psimode, phi0, verbose=verbose, plot=False, recon=recon,
			                  orbit=orbits[o])
			orbit_dat['resonance'] = resonance
			pickle.dump(orbit_dat, open(orbits[o], 'wb'))  # add resonance data to orbit_dat file
		else:
			resonance = orbit_dat['resonance']
		# r_mag = orb['OUT']['R0']  # [m]
		# Rmaj_par = np.sqrt(orb['X'] ** 2 + orb['Y'] ** 2)
		# rminor = np.sqrt((Rmaj_par - r_mag) ** 2 + orb['Z'] ** 2)
		of_l0 = (nnum * orb['CIQ']['W_PHI'] - mnum * orb['CIQ']['W_THETA']) / 2. / np.pi
		of_lm1 = (nnum * orb['CIQ']['W_PHI'] - (mnum - 1) * orb['CIQ']['W_THETA']) / 2. / np.pi
		of_lp1 = (nnum * orb['CIQ']['W_PHI'] - (mnum + 1) * orb['CIQ']['W_THETA']) / 2. / np.pi
		# orbit_r = np.mean(rminor)
		
		o_res.append(resonance['phi0'][resonance['omax']])
		# r_orb.append(orbit_r)
		f_orb[:, o] = [of_lm1, of_l0, of_lp1]
		psi_res.append(resonance['psimode'][resonance['pmax']])
		f_res.append(resonance['fepm'][resonance['fmax']])
	
	# PLOTTING
	# chosen = np.where(f_res == max(f_res))[0][0]
	chosen = closest(f_res, 75.e3)
	chosen = chosen[1]
	chosen = 0
	chosen_orb = pickle.load(open(orbits[chosen], 'rb'))
	orb = chosen_orb['orb']
	trend = chosen_orb['resonance']['trend'][:, :, chosen_orb['resonance']['omax']]
	psimode = chosen_orb['resonance']['psimode']
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
	
	# LTX cross section, ionization points, and chosen orbit
	ax1.plot(orb['OUT']['RVES'] * 100., orb['OUT']['ZVES'] * 100., 'k-')
	for p in range(len(psi_res)):
		ax1.plot(ioniz0[2, p] * 100., ioniz0[3, p] * 100., 'd', c=clrs[p % len(clrs)])  # x,y,r,z
	ax1.plot(ioniz0[2, chosen] * 100., ioniz0[3, chosen] * 100., 'rx')
	ax1.contour(orb['OUT']['X_XY'] * 100., orb['OUT']['Y_XY'] * 100., orb['OUT']['NORM_PSIXY'],
	            levels=[psimode[chosen_orb['resonance']['pmax']]], zorder=0, colors='k', linestyles='--')
	ax1.plot(orb['R'][:2000] * 100., orb['Z'][:2000] * 100., c=clrs[chosen % len(clrs)], zorder=0)
	ax1.set_xlabel('R (cm)')
	ax1.set_ylabel('Z (cm)')
	
	# spectrogram for chosen orbit
	fepm = chosen_orb['resonance']['fepm'] / 1000.
	nl = 50
	# lev=findgen(nl)/(nl-1)*(max(trend));-min(trend))+min(trend)
	ax2.contourf(psimode, fepm, trend.transpose())
	ax2.plot(psimode[chosen_orb['resonance']['pmax']], fepm[chosen_orb['resonance']['fmax']], 'rx')
	ax2.set_xlabel('r_mode (cm)')
	ax2.set_ylabel('f (kHz)')
	ax2.set_title('Power Transfer Trend (arb)')
	
	# fepm vs rmode
	for p in range(len(psi_res)):
		ax3.plot(psi_res[p], f_res[p] / 1000., 'd', c=clrs[p % len(clrs)])
	ax3.plot(psi_res[chosen], f_res[chosen] / 1000., 'rx')
	ax3.set_xlabel('mode psinorm')
	ax3.set_ylabel('fepm (kHz)')
	# power trend for chosen orbit
	ax4.plot(chosen_orb['resonance']['t'], chosen_orb['resonance']['int_power'])
	ax4.set_xlabel('time (s)')
	ax4.set_ylabel('Power Transferred (arb)')


if __name__ == '__main__':
	phase_scan_ltx()
	plt.tight_layout()
	plt.show()
	pass
