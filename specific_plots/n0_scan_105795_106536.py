import matplotlib.pyplot as plt
import numpy as np
from helpful_stuff import make_patch_spines_invisible, SimpleSignal


def n0_scan(ltx_shot):
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(12, 5))
	fig.suptitle(f'{ltx_shot}')
	fs = 9  # fontsize for legend
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
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
	ax2.yaxis.set_visible(False)
	ax2r = ax2.twinx()
	ax2r.set_ylabel('max(bdens) e12 [cm^-3]', c=clrs[0])
	ax2rr = ax2.twinx()
	ax2rr.spines["right"].set_position(("axes", 1.2))
	make_patch_spines_invisible(ax2rr)
	ax2rr.spines["right"].set_visible(True)
	ax2rr.yaxis.set_label_position('right')
	ax2rr.yaxis.set_ticks_position('right')
	ax2rr.set_ylabel('max(beam heating) [kW]', c=clrs[1])
	ax3.set_ylabel('total fraction')
	ax3.set_xlabel('n0 (ext) [#/cm^3]')
	ax4.yaxis.set_visible(False)
	ax4r = ax4.twinx()
	ax4.set_xlabel('n0 (ext) [#/cm^3]')
	ax4r.set_ylabel('Beam Heating [J]')
	
	for ax in [ax1, ax2, ax3, ax4]:
		ax.set_xscale('log')
	
	shots = [ltx_shot * 10000 + 1401 + i for i in range(15)]  # 1-15
	n0 = np.array([1.e8, 1.e9, 1.e10, 1.e11, 1.e12])
	eb_arr = np.array([10.e3, 13.e3, 16.e3])
	# n0, eb = np.meshgrid(n0_arr, eb_arr)  # same method as used in update_trdat
	ls_arr = ['-', '-.', ':']
	for i in range(len(eb_arr)):
		ax1.plot([np.nan], [np.nan], 'k', ls=ls_arr[i], label=f'{eb_arr[i] / 1000.} kV')
	ax1.legend(fontsize=fs)

	badshots = [1065361406]
	for i, (eb, ls) in enumerate(zip(eb_arr, ls_arr)):
		shotselect = shots[i * 5:i * 5 + 5]
		vbev = np.array([eb] * len(shotselect))  # [eV]
		vb = vbev / 1000.  # [kV]
		if ltx_shot == 106536:
			pb = np.array([2.87e3] * len(shotselect))
			dt = .46559 - .45834
		elif ltx_shot == 105795:
			pb = np.array([3.24e3] * len(shotselect))
			dt = .4756 - .4684  # toff-ton
		ib = pb / vbev  # [A]
		perv = ib / vbev ** 1.5
		j_tot = pb * dt  # [J] tot put into plasma
		ax1.plot(n0, perv * 1.e6, c=clrs[0])
		ax1.set_ylabel('Perveance [e-6]', c=clrs[0])
		
		ax1r.plot(n0, [p / 1000. for p in pb], c=clrs[1])
		ax1r.set_ylabel('Beam Power [kW]', color=clrs[1])
		ax1rr.plot(n0, ib, c=clrs[2])
		ax1rr.set_ylabel('Beam Current [A]', c=clrs[2])
		bden_arr = np.array([])  # [num.e12/cm^3]
		bphto_arr = np.array([])  # [kW]
		bphto_arr2 = np.array([])  # [kW]
		bpshi_arr = np.array([])  # [J]
		bplim_arr = np.array([])  # [J]
		bpte_arr = np.array([])  # [J]
		bpti_arr = np.array([])  # [J]
		cxx_arr = np.array([])  # [J]
		cxi_arr = np.array([])  # [J]
		for shot in shotselect:
			if shot in badshots:
				cxx_arr = np.append(cxx_arr, np.nan)
				cxi_arr = np.append(cxi_arr, np.nan)
				bpte_arr = np.append(bpte_arr, np.nan)
				bpti_arr = np.append(bpti_arr, np.nan)
				bpshi_arr = np.append(bpshi_arr, np.nan)
				bplim_arr = np.append(bplim_arr, np.nan)
				bphto_arr2 = np.append(bphto_arr2, np.nan)
				bphto_arr = np.append(bphto_arr, np.nan)
				bden_arr = np.append(bden_arr, np.nan)
			else:
				print(f'shot: {shot}')
				bdens = SimpleSignal(shot, '\\bdens')  # [num/cm^3]
				bphto = SimpleSignal(shot, '\\bphto')  # [W]
				bpshi = SimpleSignal(shot, '\\bpshi')  # [W]
				bplim = SimpleSignal(shot, '\\bplim')  # [W]
				bpte = SimpleSignal(shot, '\\bpte')  # [W]
				bpti = SimpleSignal(shot, '\\bpti')  # [W]
				bpcxx = SimpleSignal(shot, '\\bpcxx')  # [W]
				bpcxi = SimpleSignal(shot, '\\bpcxi')  # [W]
				
				# t_prespike = max(np.where(bpcxi.dim1 <= .473)[0])  # cut off cxi integration before spike at end of shot
				cxx_arr = np.append(cxx_arr, np.sum(bpcxx.data[1:] * (bpcxx.dim1[1:] - bpcxx.dim1[:-1])))
				cxi_arr = np.append(cxi_arr, np.sum(bpcxi.data[1:] * (bpcxi.dim1[1:] - bpcxi.dim1[: - 1])))
				bpte_arr = np.append(bpte_arr, np.sum(bpte.data[1:] * (bpte.dim1[1:] - bpte.dim1[:-1])))
				bpti_arr = np.append(bpti_arr, np.sum(bpti.data[1:] * (bpti.dim1[1:] - bpti.dim1[:-1])))
				bpshi_arr = np.append(bpshi_arr, np.sum(bpshi.data[1:] * (bpshi.dim1[1:] - bpshi.dim1[:-1])))
				bplim_arr = np.append(bplim_arr, np.sum(bplim.data[1:] * (bplim.dim1[1:] - bplim.dim1[:-1])))
				bphto_arr2 = np.append(bphto_arr2, np.sum(bphto.data[1:] * (bphto.dim1[1:] - bphto.dim1[:-1])))
				bphto_arr = np.append(bphto_arr, max(bphto.data) / 1000.)
				bden_arr = np.append(bden_arr, max(bdens.data.flatten()) / 1.e12)
		ax2r.plot(n0, bden_arr, 'o', c=clrs[0], ls=ls)
		ax2rr.plot(n0, bphto_arr, 'o', c=clrs[1], ls=ls)
		
		ax3.plot(n0, bpshi_arr / j_tot, 'o', label='shine-through', ls=ls)
		ax3.plot(n0, bplim_arr / j_tot, 'o', label='orbit loss', ls=ls)
		ax3.plot(n0, cxi_arr / j_tot, 'o', label='cx int', ls=ls)
		ax3.plot(n0, cxx_arr / j_tot, 'o', label='cx ext', ls=ls)
		
		ax4r.plot(n0, bphto_arr2, 'o', label='total', ls=ls)
		ax4r.plot(n0, bpte_arr, 'o', label='elec', ls=ls)
		ax4r.plot(n0, bpti_arr, 'o', label='ions', ls=ls)
	
	ax3.legend(fontsize=fs)
	ax4r.legend(fontsize=fs)
	plt.tight_layout()


if __name__ == '__main__':
	n0_scan(106536)
	n0_scan(105795)
	plt.show()
