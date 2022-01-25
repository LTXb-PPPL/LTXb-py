import matplotlib.pyplot as plt
import numpy as np
from transp_code.transp_classes import FBM
from helpful_stuff import SimpleSignal


def make_patch_spines_invisible(ax):
	ax.set_frame_on(True)
	ax.patch.set_visible(False)
	for sp in ax.spines.values():
		sp.set_visible(False)


# beam voltage scan (constant perveance)
# perv = 15.e-6
# vb = np.arange(10, 21) * 1000.  # 10-20kV
# ib = perv * vb ** 1.5  # A
# pb = ib * vb  # W
# vbarr = [f'{v:1.1e}' for v in vb]
# pbarr = [f'{p:1.1e}' for p in pb]
# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True)
# clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
# shots = [1036171601 + i for i in range(11)]  # P for perveance held constant
# ax1.plot(vbarr, [perv] * len(vbarr), c=clrs[0])
# ax1.set_ylabel('Perveance', c=clrs[0])
# ax1r = ax1.twinx()
# ax1r.plot(vbarr, pb, c=clrs[1])
# ax1r.set_ylabel('Beam Power [W]', color=clrs[1])
#
# bden_arr = []  # [num.e12/cm^3]
# bphto_arr = []  # [kW]
# bphto_arr2 = []  # [kW]
# bpshi_arr = []  # [J]
# bplim_arr = []  # [J]
# bpte_arr = []  # [J]
# bpti_arr = []  # [J]
# for shot in shots:
# 	print(f'shot: {shot}')
# 	bdens = SimpleSignal(shot, '\\bdens')  # [num/cm^3]
# 	bphto = SimpleSignal(shot, '\\bphto')  # [W]
# 	bpshi = SimpleSignal(shot, '\\bpshi')  # [W]
# 	bplim = SimpleSignal(shot, '\\bplim')  # [W]
# 	bpte = SimpleSignal(shot, '\\bpte')  # [W]
# 	bpti = SimpleSignal(shot, '\\bpti')  # [W]
#
# 	bpte_arr.append(np.sum(bpte.data[1:] * (bpte.dim1[1:] - bpte.dim1[:-1])))
# 	bpti_arr.append(np.sum(bpti.data[1:] * (bpti.dim1[1:] - bpti.dim1[:-1])))
# 	bpshi_arr.append(np.sum(bpshi.data[1:] * (bpshi.dim1[1:] - bpshi.dim1[:-1])))
# 	bplim_arr.append(np.sum(bplim.data[1:] * (bplim.dim1[1:] - bplim.dim1[:-1])))
# 	bphto_arr2.append(np.sum(bphto.data[1:] * (bphto.dim1[1:] - bphto.dim1[:-1])))
# 	bphto_arr.append(max(bphto.data) / 1000.)
# 	bden_arr.append(max(bdens.data.flatten()) / 1.e12)
# ax2.plot(vbarr, bden_arr, 'o-')
# ax2.set_ylabel('max(bdens) e12 [cm^-3]')
# ax2r = ax2.twinx()
# ax2r.plot(vbarr, bphto_arr, 'o-')
# ax2r.set_ylabel('max(beam heating) [kW]')
# ax3.plot(vbarr, bpshi_arr, 'o-', label='shine-through')
# ax3.plot(vbarr, bplim_arr, 'o-', label='orbit loss')
# ax4.plot(vbarr, bphto_arr2, 'o-', label='beam heating')
# ax4.plot(vbarr, bpte_arr, 'o-', label='heat to elec')
# ax4.plot(vbarr, bpti_arr, 'o-', label='heat to ions')
# ax3.set_ylabel('[J]')
# ax3.set_xlabel('Ebeam [keV]')
# ax4.set_xlabel('Ebeam [keV]')
# ax4.legend()
# ax3.legend()
#
# plt.show()
#
# a = 1

dick_notable = 0
if dick_notable:
	constant_perveance = 1
	
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(8, 5))
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	leg_fs = 10  # fontsize for legend
	fs = leg_fs + 2
	xlab = '$E_{nbi}$ [keV]'
	vbev = np.arange(10, 21) * 1000.  # 10-20kV
	xplot = vbev / 1000.  # [kV]
	shots = [1036171601 + i for i in range(11)]  # P for perveance held constant, voltage scan
	perv = np.array([15.e-6] * len(shots))
	ib = perv * vbev ** 1.5  # A
	pb = ib * vbev  # W
	
	j_tot = pb * 7.e-3  # [J] tot put into plasma
	ax1.plot(xplot, [p / 1000. for p in pb])
	ax1.set_ylabel('Beam Power [kW]', fontsize=fs)
	bden_arr = np.array([])  # [num.e12/cm^3]
	bphto_arr = np.array([])  # [kW]
	bphto_arr2 = np.array([])  # [kW]
	bpshi_arr = np.array([])  # [J]
	bplim_arr = np.array([])  # [J]
	bpte_arr = np.array([])  # [J]
	bpti_arr = np.array([])  # [J]
	cxx_arr = np.array([])  # [J]
	cxi_arr = np.array([])  # [J]
	for shot in shots:
		print(f'shot: {shot}')
		bdens = SimpleSignal(shot, '\\bdens')  # [num/cm^3]
		bphto = SimpleSignal(shot, '\\bphto')  # [W]
		bpshi = SimpleSignal(shot, '\\bpshi')  # [W]
		bplim = SimpleSignal(shot, '\\bplim')  # [W]
		bpte = SimpleSignal(shot, '\\bpte')  # [W]
		bpti = SimpleSignal(shot, '\\bpti')  # [W]
		bpcxx = SimpleSignal(shot, '\\bpcxx')  # [W]
		bpcxi = SimpleSignal(shot, '\\bpcxi')  # [W]
		
		t_prespike = max(np.where(bpcxi.dim1 <= .473)[0])  # cut off cxi integration before spike at end of shot
		cxx_arr = np.append(cxx_arr, np.sum(bpcxx.data[1:] * (bpcxx.dim1[1:] - bpcxx.dim1[:-1])))
		cxi_arr = np.append(cxi_arr,
		                    np.sum(bpcxi.data[1:t_prespike] * (
					                    bpcxi.dim1[1:t_prespike] - bpcxi.dim1[:t_prespike - 1])))
		bpte_arr = np.append(bpte_arr, np.sum(bpte.data[1:] * (bpte.dim1[1:] - bpte.dim1[:-1])))
		bpti_arr = np.append(bpti_arr, np.sum(bpti.data[1:] * (bpti.dim1[1:] - bpti.dim1[:-1])))
		bpshi_arr = np.append(bpshi_arr, np.sum(bpshi.data[1:] * (bpshi.dim1[1:] - bpshi.dim1[:-1])))
		bplim_arr = np.append(bplim_arr, np.sum(bplim.data[1:] * (bplim.dim1[1:] - bplim.dim1[:-1])))
		bphto_arr2 = np.append(bphto_arr2, np.sum(bphto.data[1:] * (bphto.dim1[1:] - bphto.dim1[:-1])))
		bphto_arr = np.append(bphto_arr, max(bphto.data) / 1000.)
		bden_arr = np.append(bden_arr, max(bdens.data.flatten()) / 1.e12)
	ax2.yaxis.set_visible(False)
	ax2r = ax2.twinx()
	ax2r.plot(xplot, bphto_arr, 'o-', c=clrs[1])
	ax2r.set_ylabel('max(Beam Heating) [kW]', fontsize=fs)
	ax2r.set_ylim(ymin=0)
	
	ax3.plot(xplot, bpshi_arr / j_tot, 'o-', label='shine-through')
	ax3.plot(xplot, bplim_arr / j_tot, 'o-', label='orbit loss')
	ax3.plot(xplot, cxi_arr / j_tot, 'o-', label='cx int')
	ax3.plot(xplot, cxx_arr / j_tot, 'o-', label='cx ext')
	ax3.set_ylabel('Total Fraction', fontsize=fs)
	ax3.set_xlabel(xlab, fontsize=fs)
	ax3.legend(fontsize=leg_fs, loc='upper left')
	
	ax4.yaxis.set_visible(False)
	ax4r = ax4.twinx()
	ax4r.plot(xplot, bphto_arr2, 'o-', label='total')
	ax4r.plot(xplot, bpte_arr, 'o-', label='elec')
	ax4r.plot(xplot, bpti_arr, 'o-', label='ions')
	ax4.set_xlabel(xlab, fontsize=fs)
	ax4r.set_ylabel('Total Beam Heating [J]', fontsize=fs)
	ax4r.legend(fontsize=leg_fs)
	
	for ax in [ax1, ax2, ax2r, ax3, ax4, ax4r]:
		ax.tick_params(labelsize=fs)
	
	plt.tight_layout()
	plt.show()

# beam scans
beam_scan = 0
if beam_scan:
	constant_power = 0
	constant_perveance = 1
	focal_length = 0
	divergence = 0
	
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(12, 5))
	# shots = 1036171401 + i for i in range(11)]  # N for neutral density scan, constant voltage, power
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	leg_fs = 10  # fontsize for legend
	fs = leg_fs + 2
	if constant_power:
		xlab = 'Ebeam [keV]'
		vbev = np.arange(10, 21) * 1000.  # 10-20kV
		xplot = vbev / 1000.  # [kV]
		shots = [1036172201 + i for i in range(11)]  # V for voltage scan, constant power
		pb = np.array([1.61e5] * len(shots))  # [W]
		ib = pb / vbev  # [A]
		perv = ib / vbev ** 1.5
	elif constant_perveance:
		xlab = 'Ebeam [keV]'
		vbev = np.arange(10, 21) * 1000.  # 10-20kV
		xplot = vbev / 1000.  # [kV]
		shots = [1036171601 + i for i in range(11)]  # P for perveance held constant, voltage scan
		perv = np.array([15.e-6] * len(shots))
		ib = perv * vbev ** 1.5  # A
		pb = ib * vbev  # W
	elif focal_length:
		xlab='Focal Length (cm)'
		foc = np.linspace(160, 200, endpoint=True, num=9)  # focal length (cm)
		shots = [1036170601 + i for i in range(9)]  # F for focal length varied
		xplot, pb, vb = foc, np.ones_like(foc)*161.e3, np.ones_like(foc)*10.7e3  # [cm], [W], [V]
		ib = pb/vb  # A
		perv = ib/vb**1.5
	elif divergence:
		xlab = 'Beam Divergence (rad)'
		div = np.linspace(.02, .1, endpoint=True, num=9)  # divergence (rad)
		shots = [1036170401 + i for i in range(9)]  # D for divergence varied
		xplot, pb, vb = div, np.ones_like(div)*161.e3, np.ones_like(div)*10.7e3  # [rad], [W], [V]
		ib = pb/vb  # A
		perv = ib/vb**1.5

	j_tot = pb * 7.e-3  # [J] tot put into plasma
	ax1.plot(xplot, perv * 1.e6, c=clrs[0])
	ax1.set_ylabel('Perveance [e-6]', c=clrs[0], fontsize=fs)
	
	ax1r = ax1.twinx()
	ax1rr = ax1.twinx()
	ax1r.spines["left"].set_position(("axes", -0.4))  # red one
	ax1rr.spines["left"].set_position(("axes", -0.2))  # green one
	make_patch_spines_invisible(ax1r)
	make_patch_spines_invisible(ax1rr)
	ax1r.spines["left"].set_visible(True)
	ax1r.yaxis.set_label_position('left')
	ax1r.yaxis.set_ticks_position('left')
	ax1rr.spines["left"].set_visible(True)
	ax1rr.yaxis.set_label_position('left')
	ax1rr.yaxis.set_ticks_position('left')
	
	ax1r.plot(xplot, [p / 1000. for p in pb], c=clrs[1])
	ax1r.set_ylabel('Beam Power [kW]', color=clrs[1], fontsize=fs)
	ax1rr.plot(xplot, ib, c=clrs[2])
	ax1rr.set_ylabel('Beam Current [A]', c=clrs[2], fontsize=fs)
	bden_arr = np.array([])  # [num.e12/cm^3]
	bphto_arr = np.array([])  # [kW]
	bphto_arr2 = np.array([])  # [kW]
	bpshi_arr = np.array([])  # [J]
	bplim_arr = np.array([])  # [J]
	bpte_arr = np.array([])  # [J]
	bpti_arr = np.array([])  # [J]
	cxx_arr = np.array([])  # [J]
	cxi_arr = np.array([])  # [J]
	for shot in shots:
		print(f'shot: {shot}')
		bdens = SimpleSignal(shot, '\\bdens')  # [num/cm^3]
		bphto = SimpleSignal(shot, '\\bphto')  # [W]
		bpshi = SimpleSignal(shot, '\\bpshi')  # [W]
		bplim = SimpleSignal(shot, '\\bplim')  # [W]
		bpte = SimpleSignal(shot, '\\bpte')  # [W]
		bpti = SimpleSignal(shot, '\\bpti')  # [W]
		bpcxx = SimpleSignal(shot, '\\bpcxx')  # [W]
		bpcxi = SimpleSignal(shot, '\\bpcxi')  # [W]
		
		t_prespike = max(np.where(bpcxi.dim1 <= .473)[0])  # cut off cxi integration before spike at end of shot
		cxx_arr = np.append(cxx_arr, np.sum(bpcxx.data[1:] * (bpcxx.dim1[1:] - bpcxx.dim1[:-1])))
		cxi_arr = np.append(cxi_arr,
		                    np.sum(bpcxi.data[1:t_prespike] * (bpcxi.dim1[1:t_prespike] - bpcxi.dim1[:t_prespike - 1])))
		bpte_arr = np.append(bpte_arr, np.sum(bpte.data[1:] * (bpte.dim1[1:] - bpte.dim1[:-1])))
		bpti_arr = np.append(bpti_arr, np.sum(bpti.data[1:] * (bpti.dim1[1:] - bpti.dim1[:-1])))
		bpshi_arr = np.append(bpshi_arr, np.sum(bpshi.data[1:] * (bpshi.dim1[1:] - bpshi.dim1[:-1])))
		bplim_arr = np.append(bplim_arr, np.sum(bplim.data[1:] * (bplim.dim1[1:] - bplim.dim1[:-1])))
		bphto_arr2 = np.append(bphto_arr2, np.sum(bphto.data[1:] * (bphto.dim1[1:] - bphto.dim1[:-1])))
		bphto_arr = np.append(bphto_arr, max(bphto.data) / 1000.)
		bden_arr = np.append(bden_arr, max(bdens.data.flatten()) / 1.e12)
	ax2.yaxis.set_visible(False)
	ax2r = ax2.twinx()
	ax2r.plot(xplot, bden_arr, 'o-', c=clrs[0])
	ax2r.set_ylabel('max(bdens) e12 [cm^-3]', c=clrs[0], fontsize=fs)
	ax2rr = ax2.twinx()
	ax2rr.spines["right"].set_position(("axes", 1.2))
	make_patch_spines_invisible(ax2rr)
	ax2rr.spines["right"].set_visible(True)
	ax2rr.yaxis.set_label_position('right')
	ax2rr.yaxis.set_ticks_position('right')
	ax2rr.plot(xplot, bphto_arr, 'o-', c=clrs[1])
	ax2rr.set_ylabel('max(beam heating) [kW]', c=clrs[1], fontsize=fs)
	ax2r.set_ylim(ymin=0)
	ax2rr.set_ylim(ymin=0)
	
	ax3.plot(xplot, bpshi_arr / j_tot, 'o-', label='shine-through')
	ax3.plot(xplot, bplim_arr / j_tot, 'o-', label='orbit loss')
	ax3.plot(xplot, cxi_arr / j_tot, 'o-', label='cx int')
	ax3.plot(xplot, cxx_arr / j_tot, 'o-', label='cx ext')
	ax3.set_ylabel('total fraction', fontsize=fs)
	ax3.set_xlabel(xlab, fontsize=fs)
	ax3.legend(fontsize=leg_fs)
	
	ax4.yaxis.set_visible(False)
	ax4r = ax4.twinx()
	ax4r.plot(xplot, bphto_arr2, 'o-', label='total')
	ax4r.plot(xplot, bpte_arr, 'o-', label='elec')
	ax4r.plot(xplot, bpti_arr, 'o-', label='ions')
	ax4.set_xlabel(xlab, fontsize=fs)
	ax4r.set_ylabel('Beam Heating [J]', fontsize=fs)
	ax4r.legend(fontsize=leg_fs)
	
	for ax in [ax1, ax1r, ax1rr, ax2, ax2r, ax2rr, ax3, ax4, ax4r]:
		ax.tick_params(labelsize=fs)

	plt.tight_layout()
	plt.show()

# ni, ne profiles avg during beam
# fig, ax = plt.subplots(figsize=(5, 4))
# shot = 1036170301
# ti = SimpleSignal(shot, '\\ti')
# te = SimpleSignal(shot, '\\te')
# tte = np.where((.463 <= te.dim2) & (.47 >= te.dim2))[0]
# tti = np.where((.463 <= ti.dim2) & (.47 >= ti.dim2))[0]
# ti_av = np.mean(ti.data[tti, :], axis=0)
# te_av = np.mean(te.data[tte, :], axis=0)
# plt.plot(te.dim1, te_av, label='te')
# plt.plot(ti.dim1, ti_av, label='ti')
# plt.legend()
# plt.xlabel('r/a')
# plt.ylabel('T (eV)')
# plt.title('<beam on>')
# plt.show()

# POH vs PBEAM
# fig, ax = plt.subplots(figsize=(5, 4))
# fig2, ax2 = plt.subplots(figsize=(5,4))
# for i, shot in enumerate([1036170301, 1034650301]):
# 	if i == 0:
# 		ls = '-'
# 	else:
# 		ax.set_prop_cycle(None)
# 		ls = '--'
# 	poh = SimpleSignal(shot, '\\poh')
# 	dvol = SimpleSignal(shot, '\\dvol')  # cm^3
# 	poh1d = np.sum(poh.data * dvol.data, axis=1)
# 	b1, = ax.plot(poh.dim2, poh1d/1000., ls=ls)
# 	pb = SimpleSignal(shot, '\\bphto')
# 	# b1, = plt.plot(poh.dim1, poh.data/1000., label='\\poh', ls=ls)
# 	b2, = ax.plot(pb.dim1, pb.data / 1000., ls=ls)
# 	if i == 0:
# 		ax.legend([b1, b2], ['\\poh', '\\bphto'])
# 	ax.set_xlabel('time (s)')
# 	ax.set_ylabel('Power (kW)')
# 	ax.set_xlim((.46, .473))
# 	ax.set_ylim((0, 400))
# 	plt.tight_layout()
#
# 	ax2.plot(pb.dim1, pb.data/poh1d, ls=ls)
# plt.show()

# STORED ENERGY
stored_energy = 1
if stored_energy:
	leg_fs = 12
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, sharex=True, figsize=(7, 5))
	desc = ['NBI and NUBEAM', 'NBI no NUBEAM', 'no NBI']
	shotset1 = [1034650301, 1034650103, 1034460103]
	shotset2 = [1036170301, 1036170103, 1036580103]
	for iss, shotset in enumerate([shotset1, shotset2]):
		for i, shot in enumerate(shotset):
			bpti = SimpleSignal(shot, '\\ti0')
			utotl = SimpleSignal(shot, '\\utotl')  # J/cm^3
			dvol = SimpleSignal(shot, '\\dvol')  # cm^3
			utot1d = np.sum(utotl.data * dvol.data, axis=1)
			if iss == 0:
				ax1.plot(bpti.dim1 * 1000., bpti.data, label=desc[i])
				ax3.plot(utotl.dim2 * 1000., utot1d)
			else:
				ax2.plot(bpti.dim1 * 1000., bpti.data)
				ax4.plot(utotl.dim2 * 1000., utot1d, label=desc[i])
	ax1.legend(loc='lower left')
	# ax2.legend(loc='lower left')
	ax3.set_xlabel('time (ms)', fontsize=leg_fs + 2)
	ax4.set_xlabel('time (ms)', fontsize=leg_fs + 2)
	ax1.set_ylabel('$T_i(0)$ (eV)', fontsize=leg_fs + 2)
	ax3.set_ylabel('$W_{tot}$ (J)', fontsize=leg_fs + 2)
	ax1.set_ylim((0, 300))
	ax2.set_ylim((0, 300))
	ax3.set_ylim((0, 800))
	ax4.set_ylim((0, 800))
	for ax in [ax1, ax2, ax3, ax4]:
		ax.set_xlim((460, 473))
		ax.axvline(463, c='k', ls='--')
		ax.axvline(470, c='k', ls='--')
		ax.tick_params(labelsize=fs)
	plt.tight_layout()
	plt.savefig('C:/Users/wcapecch/Dropbox/work_stuff/group_meetings/ti0_wtot.png')
	plt.show()

# BEAM FRACTIONS
beam_frac = 0
if beam_frac:
	fig = plt.figure(figsize=(5, 4))
	shots = [1036170301]  # , 1034650301]
	# shots = [1000020916]
	for i, shot in enumerate(shots):
		if i == 0:
			ls = '-'
		else:
			plt.gca().set_prop_cycle(None)
			ls = '--'
		bpshi = SimpleSignal(shot, '\\bpshi')
		bplim = SimpleSignal(shot, '\\bplim')
		bpcxi = SimpleSignal(shot, '\\bpcxi')
		bpcxx = SimpleSignal(shot, '\\bpcxx')
		bpth = SimpleSignal(shot, '\\bpth')
		b1, = plt.plot(bpshi.dim1, bpshi.data / 1000., label='\\bpshi', ls=ls)
		b2, = plt.plot(bplim.dim1, bplim.data / 1000., label='\\bplim', ls=ls)
		# b3, = plt.plot(bpcxi.dim1, bpcxi.data / 1000., label='\\bpcxi', ls=ls)
		b4, = plt.plot(bpcxx.dim1, bpcxx.data / 1000., label='\\bpcxx', ls=ls)
		# b5, = plt.plot(bpth.dim1, bpth.data / 1000., label='\\bpth', ls=ls)
		# plt.plot(bpshi.dim1, bpshi.data+bplim.data+bpcxi.data+bpcxx.data+bpth.data)
		if i == 0:
			plt.legend([b1, b2, b4], ['shine-through', 'orbit loss', 'cx loss (ext)'])
		# plt.legend([b1, b2, b3, b4, b5], ['\\bpshi', '\\bplim', '\\bpcxi', '\\bpcxx', '\\bpth'])
		# clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
		# tcon = [.46454, .46677, .469]
		# pbeam = 170.  # beam power [kW]
		# sh_con = np.array([.21877, .16951, .20891])*pbeam
		# pl_con = np.array([.51248, .31675, .18816])*pbeam
		# coup_con = np.array([.26876, .51374, .61292])*pbeam
		# plt.plot(tcon, sh_con, 's', c=clrs[0])
		# plt.plot(tcon, pl_con, 's', c=clrs[1])
		plt.xlabel('time (s)')
		plt.ylabel('Power (kW)')
		plt.xlim((.46, .473))
		plt.ylim((0, 60))
		plt.tight_layout()
	plt.show()

# BEAM power to ions/elec
beam_power_to_ion_elec = 0
if beam_power_to_ion_elec:
	fig = plt.figure(figsize=(5, 4))
	# shots = [1036170301, 1034650301]
	shots = [1036170301]
	for i, shot in enumerate(shots):
		if i == 0:
			ls = '-'
		else:
			plt.gca().set_prop_cycle(None)
			ls = '--'
		bpti = SimpleSignal(shot, '\\bpti')
		bpte = SimpleSignal(shot, '\\bpte')
		bphto = SimpleSignal(shot, '\\bphto')
		b3, = plt.plot(bphto.dim1, bphto.data / 1000., label='\\bphto', ls=ls)
		b1, = plt.plot(bpti.dim1, bpti.data / 1000., label='\\bpti', ls=ls)
		b2, = plt.plot(bpte.dim1, bpte.data / 1000., label='\\bpte', ls=ls)
		if i == 0:
			plt.legend([b3, b1, b2], ['beam heating tot', 'ion heating', 'elec heating'])
		plt.xlabel('time (s)')
		plt.ylabel('Power (kW)')
		plt.xlim((.46, .473))
		plt.ylim((0, 100))
		plt.axvline(.463, c='k', ls='--')
		plt.axvline(.470, c='k', ls='--')
		plt.tight_layout()
	plt.show()

# shinethrough fraction
# sh1 = SimpleSignal(1036170301, '\\sbshine_h')
# sh2 = SimpleSignal(1034650301, '\\sbshine_h')
# inj1 = SimpleSignal(1036170301, '\\sinj')
# inj2 = SimpleSignal(1034650301, '\\sinj')
# tt = np.linspace(.46, .473, num=1000)
# sh1 = np.interp(tt, sh1.dim1, sh1.data)
# sh2 = np.interp(tt, sh2.dim1, sh2.data)
# inj1 = np.interp(tt, inj1.dim1, inj1.data)
# inj2 = np.interp(tt, inj2.dim1, inj2.data)
# plt.figure(figsize=(5, 4))
# plt.plot(tt, sh1 / inj1, label='103617')
# plt.plot(tt, sh2 / inj2, label='103465')
# plt.legend()
# plt.tight_layout()
# plt.xlabel('time (s)')
# plt.ylabel('shine-through')
# plt.show()

# for sig in ['\\tauee', '\\tauea']:
# 	plt.figure()
# 	nbi1 = SimpleSignal(1036170301, sig)
# 	nbi2 = SimpleSignal(1034650301, sig)
# 	nonbi1 = SimpleSignal(1036170103, sig)
# 	nonbi2 = SimpleSignal(1034650103, sig)
# 	tt = np.linspace(.46, .473, num=1000)
# 	nbi1d = np.interp(tt, nbi1.dim1, nbi1.data)
# 	nonbi1d = np.interp(tt, nonbi1.dim1, nonbi1.data)
# 	nbi2d = np.interp(tt, nbi2.dim1, nbi2.data)
# 	nonbi2d = np.interp(tt, nonbi2.dim1, nonbi2.data)
#
# 	plt.plot(tt, nonbi1d - nbi1d, label='103617')
# 	plt.plot(tt, nonbi2d - nbi2d, label='103465')
# 	plt.xlabel('s')
# 	plt.ylabel(f'decrease in {sig} by including NUBEAM')
# 	plt.tight_layout
# 	plt.axvline(.463, c='k', ls='--')
# 	plt.axvline(.470, c='k', ls='--')
# 	plt.legend()
# plt.show()


do_sigs = 0
if do_sigs:
	# shots = [1036170301, 1034650301]
	# shots = [1036170301]#, 1036170103, 1036580103]
	# shots = [1034650301]#, 1034650103, 1034460103]
	desc = ['103617: nbi- exp and model', '103617: nbi- exp only', '103658: no nbi', '103465: nbi- exp and model',
	        '103465: nbi- exp only', '103446: no nbi']
	# shots = [1036170301, 1036170103, 1034650301, 1034650103]
	shots = [1036170103, 1036580103, 1034650103, 1034460103]
	# sigs = ['\\te', '\\te0', '\\ter_in', '\\ter_use']
	# shots = [1036170301]
	# sigs = ['\\ti0', '\\te']  # '\\bphto', '\\conde', '\\condi', '\\eheat', '\\iheat', '\\pbe', '\\pbi']
	# sigs = ['\\tauee', '\\tauea']
	# sigs = ['\\sbshine_h', '\\sinj']
	sigs = ['\\ne', '\\ni', '\\pcur', '\\pvol', '\\bdens']
	# sigs = ['\\condi','\\conde','\\CONDIWNC','\\CONDWNCE']
	# sigs = ['\\bpth', '\\bpti', '\\bpte', '\\bphto', '\\poh']
	# sigs = ['\\bphto', '\\poht']
	for sig in sigs:
		fig = plt.figure(sig, figsize=(5, 4))
		for i, sh in enumerate(shots):
			ss = SimpleSignal(sh, sig)
			if ss.ndim != 1:
				axi = fig.add_subplot(1, len(shots), i + 1, projection='3d')
				ss.plot(ax=axi)
			# plt.plot(ss.dim2, ss.data[:, 0])
			else:
				plt.plot(ss.dim1, ss.data)  # , label=desc[i])
				# plt.legend()
				plt.axvline(.463, c='k', ls='--')
				plt.axvline(.470, c='k', ls='--')
				plt.xlabel(ss.dim1units)
				plt.ylabel(sig)
				# plt.ylim((0, .007))
				plt.xlim((.46, .473))
	plt.show()

# bdens = SimpleSignal(1036170301, '\\bdens')
# vbeam = SimpleSignal(103617, '.oper_diags.ltx_nbi.source_diags.v_hvps', tree='ltx_b')
# ip = SimpleSignal(103617, '.diagnostics.magnetic.ip_rog', tree='ltx_b')
# plt.figure()
# xx, yy = np.meshgrid(bdens.dim1, bdens.dim2 * 1000.)
# plt.contourf(xx, yy, bdens.data)
# plt.rcParams.update({'font.size': 14})
# plt.colorbar(label='$n_{fi} (cm^{-3})$')
# plt.xlabel('r/a', fontsize=14)
# plt.ylabel('time (ms)', fontsize=14)
# plt.ylim((462, 480))
# plt.title('shot 103617, Ip=130kA')
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
#
# # plt.axhline(463, linestyle='--')
# # plt.axhline(471, ls='--', c='k')
# plt.show()
# a = 1
