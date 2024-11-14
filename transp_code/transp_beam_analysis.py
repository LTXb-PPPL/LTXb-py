import matplotlib.pyplot as plt
import numpy as np
from toolbox.helpful_stuff import SimpleSignal
from toolbox.helpful_stuff import smooth

clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']

# ltx_shots_nobeam = [108889, 108878, 108883, 108877, 108882, 108879, 108888, 108893, 108896, 108892, 108873]
# ltx_shots_beam = [108886, 108887, 108890, 108880, 108874, 108881, 108895, 108891, 108876, 108875, 108894]
# transp_nobeam = [1088890101]
transp_nobeam = [1088860101]
transp_beam = [1088860301]
mrks = ['', 's', 'v', 'o', 'd']
fig, (ax1, ax2) = plt.subplots(nrows=2, sharex='col', figsize=(5, 5))
hndls = []
markevery = 10
shft = 2
for i, shot in enumerate(transp_nobeam):
	c = 'k'
	ip = SimpleSignal(shot, '\\pcur')
	poh = SimpleSignal(shot, '\\poht')
	pohs = smooth(poh.data * 1.e-3, 51)
	ax1.plot(ip.dim1[shft * i:] * 1.e3, ip.data[shft * i:] * 1.e-3, c=c, marker=mrks[i], markevery=10)
	ax1.plot(poh.dim1[shft * i:] * 1.e3, pohs[shft * i:], '--', c=c, marker=mrks[i], markevery=10)
	h, = ax1.plot(np.nan, np.nan, c=c, marker=mrks[i], label=f'{shot}')
	hndls.append(h)
for i, shot in enumerate(transp_beam):
	c = 'r'
	pinj = SimpleSignal(shot, '\\pinj')
	pb = SimpleSignal(shot, '\\bpcap')
	pbe = SimpleSignal(shot, '\\pbe')  # Watts/cm^3
	dvol = SimpleSignal(shot, '\\dvol')  # cm^3
	pbe1d = np.sum(pbe.data * dvol.data, axis=1)
	ip = SimpleSignal(shot, '\\pcur')
	poh = SimpleSignal(shot, '\\poht')
	pohs = smooth(poh.data * 1.e-3, 51)
	ax2.plot(pinj.dim1[shft * i:] * 1.e3, pinj.data[shft * i:] * 1.e-3, c=c, marker=mrks[i], markevery=10)
	ax2.plot(pb.dim1[shft * i:] * 1.e3, pb.data[shft * i:] * 1.e-3, '--', c=c, marker=mrks[i], markevery=10)
	ax2.plot(pbe.dim2[shft * i:] * 1.e3, pbe1d[shft * i:] * 1.e-3, '-.', c=c, marker=mrks[i], markevery=10)
	ax1.plot(ip.dim1[shft * i:] * 1.e3, ip.data[shft * i:] * 1.e-3, c=c, marker=mrks[i], markevery=10)
	ax1.plot(poh.dim1[shft * i:] * 1.e3, pohs[shft * i:], '--', c=c, marker=mrks[i], markevery=10)
	h, = ax1.plot(np.nan, np.nan, c=c, marker=mrks[i], label=f'{shot}')
	hndls.append(h)

fs, fs2 = 12, 10
ax2.set_ylabel('$P_{NBI}$ (kW)', fontsize=fs)
ax2.set_xlabel('time (ms)', fontsize=fs)
# ax1.set_ylabel('$I_p$ (kA) -\n$P_{OH}$ (kW)')
ax2.set_xlim((450, 480))
for ax in [ax1, ax2]:
	ax.tick_params(labelsize=fs2)
l1, = ax1.plot(np.nan, np.nan, 'k-', label='$I_p$ (kA)')
l2, = ax1.plot(np.nan, np.nan, 'k--', label='$P_{OH} (kW)$')
leg1 = ax1.legend(handles=[l1, l2], loc='upper center', labelcolor='linecolor', frameon=False, fontsize=fs2)
l3, = ax1.plot(np.nan, np.nan, 'k', label='No Beam', ls='none')
l4, = ax1.plot(np.nan, np.nan, 'r', label='With Beam', ls='none')
leg2 = ax1.legend(handles=[l3, l4], loc='upper left', labelcolor='linecolor', frameon=False, fontsize=fs2,
                  handlelength=0)
ax1.legend(handles=hndls, loc='upper right', frameon=False, fontsize=fs2)
ax1.add_artist(leg1)
ax1.add_artist(leg2)
l5, = ax2.plot(np.nan, np.nan, 'k-', label='injected')
l6, = ax2.plot(np.nan, np.nan, 'k--', label='captured')
l7, = ax2.plot(np.nan, np.nan, 'k-.', label='to electrons')
ax2.legend(handles=[l5, l6, l7], fontsize=fs2, frameon=False)

plt.tight_layout()
plt.show()

# shts = [1083040301, 1084350301, 1083690301]
# # sigs = ['\\pinj', '\\bpcap', '\\bplim', '\\bdens', '\\bphto', '\\pbe', '\\bpshi']
# psigs = ['\\pinj', '\\bpcap', '\\bplim', '\\bphto', '\\bpshi']
#
# fig, ax = plt.subplots()
# for i, shot in enumerate(shts):
# 	c = clrs[i]
# 	pinj = SimpleSignal(shot, '\\pinj')
# 	pb = SimpleSignal(shot, '\\bpcap')
# 	pbe = SimpleSignal(shot, '\\pbe')  # Watts/cm^3
# 	dvol = SimpleSignal(shot, '\\dvol')  # cm^3
# 	pbe1d = np.sum(pbe.data * dvol.data, axis=1)
# 	if i == 0:
# 		ax.plot(pinj.dim1 * 1.e3, pinj.data * 1.e-3, c=c, label='injected')
# 		ax.plot(pb.dim1 * 1.e3, pb.data * 1.e-3, '--', c=c, label='captured')
# 		ax.plot(pbe.dim2 * 1.e3, pbe1d * 1.e-3, '-.', c=c, label='to electrons')
# 	else:
# 		ax.plot(pinj.dim1 * 1.e3, pinj.data * 1.e-3, c=c)
# 		ax.plot(pb.dim1 * 1.e3, pb.data * 1.e-3, '--', c=c)
# 		ax.plot(pbe.dim2 * 1.e3, pbe1d * 1.e-3, '-.', c=c)
# 	poh = SimpleSignal(shot, '\\poht')
# 	ax1.plot(ip.dim1 * 1.e3, ip.data * 1.e-3, c=c)
# 	pohs = smooth(poh.data * 1.e-3, 51)
# 	ax1.plot(poh.dim1 * 1.e3, pohs, '--', c=c)
#
# ax.set_xlim((460,480))
# ax.set_xlabel('time (ms)')
# ax.set_ylabel('Power (kW)')
# plt.tight_layout()
# plt.show()

if __name__ == '__main__':
	pass
