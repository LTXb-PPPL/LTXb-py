import matplotlib.pyplot as plt
import numpy as np
from helpful_stuff import SimpleSignal


def plot_sig(sig='.inputs.ter', shot=1000030301, ax=None):
	dat = SimpleSignal(shot, sig, tree='transp_ltx')
	dat.plot(ax=ax, label=f'{sig}, shot: {shot}')
	# t = .4583
	# ii = max(np.where(dat.dim2 <= t)[0])
	# plt.plot(dat.dim1, dat.data[ii, :], label=shot)


if __name__ == '__main__':
	shts = [1065361802, 1057951802]
	for i, sh in enumerate(shts):
		fig, axs = plt.subplots(nrows=2, ncols=3, sharex='col')
		fig.suptitle(f'shot: {sh}')
		for ii, sig in enumerate(['\\ashaf', '\\raxis', '\\pcur', '\\pvol', '\\bplim', '\\ne']):
			dat = SimpleSignal(sh, sig, tree='transp_ltx')
			if sig == '\\ne':
				axs.flatten()[ii].plot(dat.dim2, dat.data[:, 0])  # plot core ne
			else:
				axs.flatten()[ii].plot(dat.dim1, dat.data)
			axs.flatten()[ii].set_ylabel(f'{sig}')
			if i == 0:
				axs.flatten()[ii].axvline(.458, ls='--',c='k')
				axs.flatten()[ii].axvline(.466, ls='--', c='k')
			else:
				axs.flatten()[ii].axvline(.468, ls='--',c='k')
				axs.flatten()[ii].axvline(.476, ls='--', c='k')
	plt.tight_layout()
	plt.show()
