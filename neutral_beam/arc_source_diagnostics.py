import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from toolbox.helpful_stuff import get_tree_conn, get_data

"""
Diagnosing strange arc source behavior
Data from 5/11/22
"""


def arcs_11may22():
	csvdir = 'C:/Users/wcapecch/Dropbox/work_stuff/nbi/arc_source_testing_data/'
	shots1 = 508615 + np.array([23, 25, 26, 30, 32])  # Iarc=400
	shots2 = 508615 + np.array([45, 48, 50, 52, 53])  # Iarc=500
	shots1_fn = [5, 6, 7, 9, 10]
	shots2_fn = [12, 13, 14, 15, 16]
	
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex='col')
	ax1r, ax3r = ax1.twinx(), ax3.twinx()
	ax3.set_xlabel('time (s)')
	ax4.set_xlabel('time (s)')
	ax1.set_ylabel('Iarc (A)')
	ax1r.set_ylabel('Varc (V)')
	ax2.set_ylabel('Iarcps (V)')
	ax3.set_ylabel('Iarc (A)')
	ax3r.set_ylabel('Varc (V)')
	ax4.set_ylabel('Iarcps (V)')
	for ax in [ax2, ax4]:
		ax.set_xlim((-.003, .008))
	
	for sh, ifn in zip(shots1, shots1_fn):
		print(sh)
		t = get_tree_conn(sh, treename='ltx_nbi')
		(tiarc, iarc) = get_data(t, '.source_diags.i_arc')
		(tarc, varc) = get_data(t, '.source_diags.v_arc')
		iarc = np.interp(tarc, tiarc, iarc)
		ax1.plot(tarc, iarc)
		ax1r.plot(tarc, varc)
		csvfn = f'{csvdir}NewFile{ifn}.csv'
		csv = pd.read_csv(csvfn, header=1, index_col=False)
		ax2.plot(csv['Second'], csv['Volt'])
	for sh, ifn in zip(shots2, shots2_fn):
		t = get_tree_conn(sh, treename='ltx_nbi')
		(tiarc, iarc) = get_data(t, '.source_diags.i_arc')
		(tarc, varc) = get_data(t, '.source_diags.v_arc')
		iarc = np.interp(tarc, tiarc, iarc)
		ax3.plot(tarc, iarc)
		ax3r.plot(tarc, varc)
		csvfn = f'{csvdir}NewFile{ifn}.csv'
		csv = pd.read_csv(csvfn, header=1, index_col=False)
		ax4.plot(csv['Second'], csv['Volt'])
	plt.tight_layout()
	plt.show()


def flat_iarc():
	# data taken 12May22
	shots = 508723 + np.array([18, 21, 22, 24, 26, 29, 34])
	newfilenums = [0, 1, 2, 3, 4, 5, 6]
	
	fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(nrows=2, ncols=2, figsize=(10, 5))
	ax21r, ax22r = ax21.twinx(), ax22.twinx()
	# for sh, nfn in zip(shots, newfilenums):
	for sh in shots:
		print(sh)
		t = get_tree_conn(sh, treename='ltx_nbi')
		(tiarc, iarc) = get_data(t, '.source_diags.i_arc')
		(tarc, varc) = get_data(t, '.source_diags.v_arc')
		iarc = np.interp(tarc, tiarc, iarc)
		ax11.plot(tarc, varc, label=sh)
		ax12.plot(tarc, iarc)
	for sh in 508723 + np.array([58, 60]):
		print(sh)
		t = get_tree_conn(sh, treename='ltx_nbi')
		(tiarc, iarc) = get_data(t, '.source_diags.i_arc')
		(tarc, varc) = get_data(t, '.source_diags.v_arc')
		iarc = np.interp(tarc, tiarc, iarc)
		(tihv, ihv) = get_data(t, '.source_diags.i_hvps')
		(tbeam, vbeam) = get_data(t, '.source_diags.v_hvps')
		ibeam = np.interp(tbeam, tihv, ihv)
		ax21.plot(tarc, varc, label=sh)
		ax21r.plot(tbeam, vbeam, '--')
		ax22.plot(tarc, iarc)
		ax22r.plot(tbeam, ibeam, '--')
	
	ax11.legend()
	ax21.legend()
	for ax in [ax21, ax22]:
		ax.set_xlabel('time (s)')
	for ax in [ax11, ax21]:
		ax.set_ylabel('$V_{arc}$ (V)')
	for ax in [ax12, ax22]:
		ax.set_ylabel('$I_{arc}$ (A)')
	ax21r.set_ylabel('$V_{beam}$ (V)')
	ax22r.set_ylabel('$I_{beam}$ (A)')
	for aa in [ax11, ax12, ax21, ax22]:
		aa.grid(b=True, which='both')
	plt.tight_layout()
	plt.show()


if __name__ == '__main__':
	# arcs_11may22()
	flat_iarc()
