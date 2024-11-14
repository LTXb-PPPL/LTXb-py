import matplotlib.pyplot as plt
import numpy as np
from toolbox.ionization_cross_sections import i_ch_ex_data, iimpact_data, tabulated2_eimpact_ionization

clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']


def li_cx(show_fns=False):
	'''This routine plots the cross sections taken from ADAS of H^0+Li^3+ to various states of H^+ + Li^2+'''
	reac = '$Li^{3+} + H^0 -> Li^{2+}+H^{-1}$'
	direc = 'Z:/PycharmProjects/LTXb-py/neutral_beam/li_beam_fueling/'
	fns = ['ext#h0_arf07#li3', 'qcx#h0_en2_kvi#li3', 'qcx#h0_gyt#li3']
	for fn in fns:
		en, cx, en_units, cx_units, nmin, nmax = get_total_cx_adas(f'{direc}{fn}.dat')
		if show_fns:
			lbl = f'{reac} : n={nmin}-{nmax} : {fn}'
		else:
			lbl = f'{reac} : n={nmin}-{nmax}'
		plt.plot(en, cx, label=lbl)
	hcx_en, hcx_xs = i_ch_ex_data()
	hii_en, hii_xs = iimpact_data()
	hei_en, hei_xs = tabulated2_eimpact_ionization(200.)
	hcx_en, hii_en, hei_en = hcx_en / 1000., hii_en / 1000., hei_en / 1000.  # convert eV to keV/amu
	plt.plot(hcx_en, hcx_xs, '--', label='H-H cx')
	plt.plot(hii_en, hii_xs, '--', label='H-H i-impact')
	plt.plot(hei_en, hei_xs, '--', label='H-H e-impact (Te=200eV)')
	plt.xlabel(f'energy {en_units}')
	plt.ylabel(f'$\sigma$ {cx_units}')
	plt.xscale('log')
	plt.xlim(1.e-3, 1.e3)
	plt.legend()


def get_total_cx_adas(fn):
	'''using knowledge of structure here- 9 columns, energy and total cx arrays first followed by partials'''
	print(f'\n{fn}')
	with open(fn, 'r') as file:
		en, cx = [], []
		lines = file.readlines()
		lnmin, lnmax = [i for i in range(len(lines)) if 'nmin' in lines[i]], [i for i in range(len(lines)) if
		                                                                      'nmax' in lines[i]]
		nmin, nmax = int(lines[lnmin[0]].split()[0]), int(lines[lnmax[0]].split()[0])
		inum_en = [i for i in range(len(lines)) if 'number of energies' in lines[i] and not lines[i].startswith('C')]
		ien = [i for i in range(len(lines)) if 'keV/amu' in lines[i] and not lines[i].startswith('C')]
		icx = [i for i in range(len(lines)) if 'total xsects' in lines[i] and not lines[i].startswith('C')]
		for j in range(len(ien)):
			nen = int(lines[inum_en[j]].split()[0])
			lines[ien[j]] = lines[ien[j]].replace('D', 'E')  # sometimes data reported as 1D-15 instead of 1E-15
			lines[icx[j]] = lines[icx[j]].replace('D', 'E')
			en.extend([float(e) for e in lines[ien[j]].split()[:nen]])
			cx.extend([float(x) for x in lines[icx[j]].split()[:nen]])
		en_units = lines[ien[j]].split()[-1]
		cx_units = lines[icx[j]].split()[-1]
	return en, cx, en_units, cx_units, nmin, nmax


if __name__ == '__main__':
	li_cx()
	plt.show()
