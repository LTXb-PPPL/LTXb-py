import matplotlib.pyplot as plt
import numpy as np
from helpful_stuff import SimpleSignal


def plot_sig(sig='.inputs.ter', shot=1000030301, ax=None):
	dat = SimpleSignal(shot, sig, tree='transp_ltx')
	# dat.plot(ax=ax, label=f'{sig}, shot: {shot}')
	t = .4583
	ii = max(np.where(dat.dim2 <= t)[0])
	plt.plot(dat.dim1, dat.data[ii, :], label=shot)


if __name__ == '__main__':
	shts = [1065360101, 1059520301]
	for sh in shts:
		plot_sig(sig='.inputs.qpr', shot=sh)
	plt.legend()
	plt.show()
