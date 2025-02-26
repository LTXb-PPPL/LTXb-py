import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import os


def get_home_direc():
	pppl = 'Z:/PycharmProjects/LTXb-py/'
	home = 'C:/Users/willi/PycharmProjects/LTXb-py/'
	if os.path.exists(pppl):
		return pppl
	elif os.path.exists(home):
		return home
	else:
		print('no home direc found')
		return None


def smooth_onto(xonto, xsig, sig, xinterval):
	avsig = np.zeros_like(xonto)
	for i in range(len(avsig)):  # go through each pt in new array
		ii = np.where((xsig <= xonto[i] + xinterval / 2.) & (xsig > xonto[i] - xinterval / 2.))
		if len(ii[0]) > 0:
			avsig[i] = np.nanmean(sig[ii])
	return avsig


clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
matplotlib.use('TkAgg')  # allows plotting in debug mode

hdirec = get_home_direc()
fps = [f'{hdirec}npa/calibration/huzzah.dat', f'{hdirec}npa/calibration/huzzah2.dat']

fig, [ax1, ax2, ax3] = plt.subplots(nrows=3, sharex='col')
dig = 1000.  # Hz digitization rate
t1 = 4  # seconds of files to use to compute offsets
vmax = 12500  # max V for plotting
vstep = 50  # V window to average data
phv = np.arange(0, vmax + vstep, vstep)

for ifp, fp in enumerate(fps):
	with open(fp, 'r') as file:
		lines = file.readlines()
		dat = np.array([[float(a) for a in l.split('\t')] for l in lines])
	
	if len(dat[0, :]) > len(dat[:, 0]):
		dat = dat.transpose()
	t = np.arange(len(dat[:, 0])) / dig
	j = np.where(t < t1)
	offs = [np.mean(dat[:, i][j]) for i in range(8)]
	for i in range(7):  # flip channeltron signals and remove offset
		dat[:, i] = -dat[:, i] + offs[i]
	dat[:, 7] -= offs[7]  # remove ionsrc offset
	hv = dat[:, 7] * 3000.  # scale by ross divider ratio to get HV value on ion source
	ihv = np.argsort(hv)  # indexes of increasing hv
	hv = hv[ihv]  # resort hv array
	smdat = np.zeros((len(phv), 7))
	for i in range(7):
		dat[:, i] = dat[ihv, i]  # sort rows so interp works
		smdat[:, i] = smooth_onto(phv, hv, dat[:, i], vstep)  # smooth signals
	
	if ifp == 0:
		for i in range(7):
			ax1.plot(hv, dat[:, i], c=clrs[i], label=f'sig{i} raw')
		sigs = smdat
	else:
		for i in range(7):
			ax2.plot(hv, dat[:, i], c=clrs[i], label=f'sig{i} raw')
		sigs2 = smdat
		for i in range(7):
			ax3.plot(phv, np.average([sigs[:, i], sigs2[:, i]], axis=0), c=clrs[i])

ax3.set_ylabel('av sig (V)')
ax3.set_xlabel('HV-ion source (kV)')
ax1.set_ylabel('raw & interp sig (V)')
ax2.set_ylabel('raw & interp sig (V)')
plt.show()

if __name__ == '__main__':
	pass
