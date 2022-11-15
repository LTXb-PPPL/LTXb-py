import matplotlib.pyplot as plt
import numpy as np

from toolbox.helpful_stuff import SimpleSignal

shot = 1065362000  # 'T' shot group
rtan = np.array([19., 28., 36., 44., 53.])
''' Orientations:
0 = normal
1 = rev-Ip
2 = rev-Bt
3 = rev-Ip&Bt
'''
orients = [0, 1, 2, 3]
desc = ['normal', 'rev-Ip', 'rev-Bt', 'rev-Ip&Bt']
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
fs = 12
fig, (ax, ax2) = plt.subplots(nrows=2, sharex='col')
ax.set_ylabel('max(bdens) e12 [cm^-3]', fontsize=fs)
ax2.set_ylabel('Beam Heating [J]', fontsize=fs)

for o in orients:
	shots = shot + np.array(np.arange(5) + (len(rtan) * o + 1))
	maxbdens, maxheat = [], []
	for sh in shots:
		try:
			bdens = SimpleSignal(sh, '\\bdens')  # [num/cm^3]
			maxbdens.append(max(bdens.data.flatten()) / 1.e12)
			bphto = SimpleSignal(sh, '\\bphto')  # [W]
			maxheat.append(np.sum(bphto.data[1:] * (bphto.dim1[1:] - bphto.dim1[:-1])))
		except AttributeError:
			print(f'trouble fetching data for shot {sh}')
			maxbdens.append(np.nan)
			maxheat.append(np.nan)
	ax.plot(rtan, maxbdens, 'o-', c=clrs[o], label=desc[o])
	ax2.plot(rtan, maxheat, 'o-', c=clrs[o])

ax.legend()
plt.show()

if __name__ == '__main__':
	pass
