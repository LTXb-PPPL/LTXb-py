import matplotlib.pyplot as plt
import numpy as np
from toolbox.helpful_stuff import ltx_limiter, read_eqdsk3, get_tree_conn, get_data, smooth, make_patch_spines_invisible


def shot106536_plasma_params():
	eq = 'Z:/transp/t106536/106536R02_05.eqdsk'
	deq = read_eqdsk3(eq)
	
	lr, lz, _, _ = ltx_limiter()
	
	fig2, (ax2, ax1) = plt.subplots(ncols=2, figsize=(10, 4))
	ax1.contour(deq['x_xy'], deq['y_xy'], deq['psixy'], levels=np.linspace(deq['psimag'], deq['psilim'], num=7),
	            linestyles='--', colors='k')
	ax1.plot(lr, lz, 'k-')
	ax1.plot(deq['rlimit'], deq['zlimit'], 'r-')
	ax1.axis('equal')
	ax1.set_ylim((-1.1 * max(lz), 1.1 * max(lz)))
	ax1.axis('off')
	
	prefix = '.oper_diags.ltx_nbi'
	t = get_tree_conn(106536, treename='ltx_b')
	(tip, ip) = get_data(t, '.diagnostics.magnetic.ip_rog')
	(tne, ne) = get_data(t, '.diagnostics.microwave.interferom.phase_comp.ne_l')  # [m^-2]
	(tibeam, ibeam) = get_data(t, f'{prefix}.source_diags.i_hvps')
	(tbeam, vbeam) = get_data(t, f'{prefix}.source_diags.v_hvps')
	# (tacc, iacc) = get_data(t, f'{prefix}.source_diags.i_accel_grid')
	# (tdec, idec) = get_data(t, f'{prefix}.source_diags.i_decel_grid')
	ibeam = np.interp(tbeam, tibeam, ibeam)
	pbeam = ibeam * vbeam / 1000.  # kW
	
	ne = smooth(ne)
	
	clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
	# ax2.yaxis.set_visible(False)
	ax2r = ax2.twinx()
	ax2rr = ax2.twinx()
	ax2rr.spines["right"].set_position(("axes", 1.2))
	make_patch_spines_invisible(ax2rr)
	ax2rr.spines["right"].set_visible(True)
	ax2rr.yaxis.set_label_position('right')
	ax2rr.yaxis.set_ticks_position('right')
	
	fs = 16
	ax2.set_ylabel('$I_p$ (kA)', c=clrs[0], fontsize=fs)
	ax2r.spines['left'].set_color(clrs[0])
	ax2.tick_params(axis='y', colors=clrs[0])
	ax2r.set_ylabel('$\int{n_e dl}$ (e19 $m^{-2}$)', c=clrs[1], fontsize=fs)
	ax2r.spines['right'].set_color(clrs[1])
	ax2r.tick_params(axis='y', colors=clrs[1])
	ax2rr.set_ylabel('$P_{NBI}$ (kW)', c=clrs[2], fontsize=fs)
	ax2rr.spines['right'].set_color(clrs[2])
	ax2rr.tick_params(axis='y', colors=clrs[2])
	ax2.set_xlabel('time (ms)', fontsize=fs)
	ax2.plot(tip * 1000, -ip / 1000., c=clrs[0])
	ax2r.plot(tne * 1000, ne * 1.e-19, c=clrs[1])
	ax2rr.plot(tbeam * 1000, pbeam, c=clrs[2])
	ax2.axvline(464, linestyle='--', c='k')
	ax2.set_xlim((455, 470))
	for ax in [ax2, ax2r, ax2rr]:
		ax.set_ylim(bottom=0)
		ax.tick_params(labelsize=fs-2)
	plt.tight_layout()


if __name__ == '__main__':
	shot106536_plasma_params()
	plt.show()
