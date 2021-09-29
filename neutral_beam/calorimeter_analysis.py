import matplotlib.pyplot as plt
import numpy as np
from helpful_stuff import SimpleSignal
import lvm_read
from bills_LTX_MDSplus_toolbox import *
import pickle


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
# norm_cal = np.zeros((len(shots), 5))  # index of [shot, cal]
cal_arr = np.zeros((len(shots), 5))
ja = []  # shot arrays for joules

dbug = 0
if dbug:
	fig2, dax = plt.subplots()

for ish, shot in enumerate(shots):
	lvm_file = f'{dir}{shot}.lvm'
	if os.path.exists(lvm_file):
		print(shot)
		if shot > 200000:
			tree = 'ltx_nbi'
			nbi_only = True
			prefix = ''
		else:
			tree = 'ltx_b'
			nbi_only = False
			prefix = '.oper_diags'
			
		tree = get_tree_conn(shot, treename=tree)
		(ti, ib) = get_data(tree, f'{prefix}.ltx_nbi.source_diags.i_hvps')
		(tv, vb) = get_data(tree, f'{prefix}.ltx_nbi.source_diags.v_hvps')
		ib = np.interp(tv, ti, ib)  # get ib onto vb axis
		t_beamon = np.where(vb > 5000)  # only look at where beam is above 5kV
		if len(t_beamon[0]) < 10:
			ja.append(np.nan)
			cal_arr[ish, :] = np.nan
		else:
			pb = ib * vb  # beam power [W]
			tot_joules = np.sum([pb[i] * (tv[i] - tv[i - 1]) for i in np.arange(1, len(tv))])
			ja.append(tot_joules)
			
			lvm = lvm_read.read(lvm_file, dump_file=False, read_from_pickle=False)
			lvm0 = lvm[0]
			rtd_nam = lvm0['Channel names'][1:5]  # u1-4,l1-4,d1-11
			for i in np.arange(1, 5):
				cal_arr[ish, i - 1] = (max(lvm0['data'][:, i]) - min(lvm0['data'][:, i]))
	else:
		ja.append(np.nan)
		cal_arr[ish, :] = np.nan

c_ps = 385.  # J/kg/C
m_copp = 1.5834  # kg
dt_predict = [j / c_ps / m_copp for j in ja]
fig, ax1 = plt.subplots(sharey=True)#, figsize=(12, 6))
ax1.plot(ja, cal_arr[:, 0:4] * 1.e3, 'o', label=rtd_nam[0:4])  # calorimeter 0-4
ax1.set_ylabel('$\Delta T$ ($\degree$C)')
ax1.set_title('Calorimeter')
# for ax in [ax1, ax2, ax3]:
# 	ax.set_xlim(left=7.5)
# 	ax.set_ylim(top=2.)

fig2, (ax11, ax22, ax33) = plt.subplots(ncols=3, sharey=True, figsize=(12, 6))
ax11.plot(ja, cal_arr[:, 0:4], 'o', label=rtd_nam[0:4])  # upper scrapers
# ax22.plot(ja, cal_arr[:, 4:8], 'o', label=rtd_nam[4:8])  # lower scrapers
# ax33.plot(ja, cal_arr[:, 8:19], 'o', label=rtd_nam[8:19])  # dump
ax11.set_ylabel('$\Delta T$ ($\degree$C)')
ax22.set_xlabel('$E_{tot}$')
ax11.set_title('Upper')
ax22.set_title('Lower')
ax33.set_title('Dump')
for ax in [ax11, ax22, ax33]:
	pass

# ax1t = ax1.twiny()
# ax1t.plot(pa, np.transpose(np.transpose(np.outer(np.ones(3), rtd_arr[2, 0:4])) * ja / ja[2]), '--')
# ax1t.set_xlabel('Perv (e-6)')
# ax2t = ax2.twiny()
# ax2t.plot(pa, np.transpose(np.transpose(np.outer(np.ones(3), rtd_arr[2, 4:8])) * ja / ja[2]), '--')
# ax2t.set_xlabel('Perv (e-6)')
# ax3t = ax3.twiny()
# ax3t.plot(pa, np.transpose(np.transpose(np.outer(np.ones(3), rtd_arr[2, 8:19])) * ja / ja[2]), '--')
# ax3t.set_xlabel('Perv (e-6)')

ax1.legend(fontsize=10)
ax2.legend(fontsize=10)
ax3.legend(ncol=2, fontsize=10)
ax11.legend(fontsize=10)
ax22.legend(fontsize=10)
ax33.legend(ncol=2, fontsize=10)
plt.show()
