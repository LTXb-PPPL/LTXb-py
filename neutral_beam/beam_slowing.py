import matplotlib.pyplot as plt
import numpy as np

# constants
me = 9.10938e-31  # electron mass [kg]
mp = 1.67262e-27  # proton mass [kg]
qe = 1.60218e-19  # electron charge [C]
eps0 = 8.85419e-12  # vacuum permittivity [A^2*s^4/m^3/kg]


def plot_fractional_energy_transfer(eb=20000, mb=1., mi=1., interps=None):
	t_crit = .342 * (me / (mi * mp)) ** (1 / 3) / (1 + (mb / mi)) ** (2 / 3) * eb  # [eV]
	print(f'T_crit: {t_crit}eV')
	tt = np.linspace(0, 2 * t_crit)
	fe = t_crit / (t_crit + tt)  # fraction of energy to electrons
	fi = tt / (t_crit + tt)  # fraction to ions
	if interps is not None:
		f_e, f_i = np.interp(interps, tt, fe), np.interp(interps, tt, fi)
		print(f'Fractions at {interps} are to elec: {f_e}, ions: {f_i}')
	# plt.axvline(interps, c='k', ls='--')
	plt.plot(tt, fe, label='fraction to e')
	plt.plot(tt, fi, label='fraction to i')
	plt.legend()
	plt.xlabel('T (eV)')


# plt.show()


def plot_collision_freqs(axs=None, ls='-', showlegend=True):
	if axs is None:
		fig, axs = plt.subplots(2, 2, sharex='col')
		y0max = y1max = 0
	else:
		axs[0, 0].set_prop_cycle(None)
		axs[1, 0].set_prop_cycle(None)
		axs[0, 1].set_prop_cycle(None)
		axs[1, 1].set_prop_cycle(None)
		y0max = axs[0, 0].get_ylim()[1]
		y1max = axs[1, 0].get_ylim()[1]
	axr = axs[1, 0].twinx()
	mr = Ab * mp * mp / (Ab * mp + mp)  # reduced mass [kg]
	
	def get_nu_be(vb_func):
		vte = np.sqrt(2. * te * qe / me)  # electron thermal velocity [m/s]
		nu_be_func = Zb ** 2 * qe ** 4 * ne * lnlamb / (
				4 * np.pi * eps0 ** 2 * me * Ab * mp * (vb_func ** 3 + 1.3 * vte ** 3))  # Friedberg eq 9.63
		return nu_be_func
	
	def get_nu_bi(vb_func):
		vti = np.sqrt(2. * ti * qe / mp)  # ion thermal velocity [m/s]
		nu_bi_f2 = Zb ** 2 * qe ** 4 * ni * lnlamb / (
				4 * np.pi * eps0 ** 2 * mr * Ab * mp * (vb_func ** 3 + 1.3 * vti ** 3))  # Friedberg eq 9.65
		return nu_bi_f2
	
	eb_arr = np.linspace(0, eb * 2., endpoint=True, num=5000)
	eb_arr = eb_arr[1:]  # drop E=0
	vb_arr = np.sqrt(2. * eb_arr * qe / (Ab * mp))
	vb = np.sqrt(2. * eb * qe / (Ab * mp))  # beam ion velocity [m/s]
	nu_be = get_nu_be(vb_arr)
	nu_bi = get_nu_bi(vb_arr)
	nu_be_in = get_nu_be(vb)
	nu_bi_in = get_nu_bi(vb)
	
	ecrit = Ab * mp / 2. * ((1.3 * ((2. * te * qe / me) ** 1.5 * ni / mr - (2. * ti * qe / mp) ** 1.5 * ne / me)) / (
			ne / me - ni / mr)) ** (2. / 3)
	ecrit_ev = ecrit / qe
	vcrit = np.sqrt(2. * ecrit / (Ab * mp))  # critical velocity [m/s]
	w0 = vb / vcrit
	t_therm = np.log(1. + w0 ** 3) / 3. / nu_be_in  # thermalization time see Friedberg eq. 9.71
	print('\nTe: {}eV  Ti: {}eV'.format(te, ti))
	print('nu_be @ {}keV: {:.2e}'.format(eb / 1000., nu_be_in))
	print('nu_bi @ {}keV: {:.2e}'.format(eb / 1000., nu_bi_in))
	print('Ecrit (keV) = {:.2f}'.format(ecrit_ev / 1000.))
	print('t_therm: {}s'.format(t_therm))
	
	pmax = max([nu_be_in, nu_bi_in, y0max])
	nube, = axs[0, 0].plot(eb_arr / 1000., nu_be, ls=ls, label='nu_be')
	nubi, = axs[0, 0].plot(eb_arr / 1000., nu_bi, ls=ls, label='nu_bi')
	axs[0, 0].axvline(eb / 1000., c='k', alpha=0.5, ls='--')
	axs[0, 0].annotate('Ebeam', (eb / 1000., 0.5 * pmax))
	axs[0, 0].axvline(ecrit_ev / 1000., c='k', alpha=0.5, ls='-.')
	axs[0, 0].annotate('Ecrit', (ecrit_ev / 1000., 0.5 * pmax))
	axs[0, 0].set_ylabel('Collision Freq [s^-1]')
	axs[0, 0].set_ylim((0, 1.3 * pmax))
	if showlegend:
		axs[0, 0].legend([nube, nubi], ['nu_be', 'nu_bi'], fontsize=9)
	
	dpedt = nu_be * Ab * mp * vb_arr  # momentum transfer to electrons (see Friedberg eq 9.26)
	dpidt = nu_bi * Ab * mp * vb_arr  # [kgm/s/s]
	dpdt = dpedt + dpidt  # total momentum transfer rate [kgm/s/s]
	pmax = max([nu_be_in * Ab * mp * vb, nu_bi_in * Ab * mp * vb, y1max])
	e10, = axs[1, 0].plot(eb_arr / 1000., dpedt, ls=ls)
	i10, = axs[1, 0].plot(eb_arr / 1000., dpidt, ls=ls)
	t10, = axs[1, 0].plot(eb_arr / 1000., dpdt)
	axs[1, 0].set_xlabel('Beam Energy [keV]')
	axs[1, 0].set_ylabel('dp/dt')
	axs[1, 0].set_ylim((0, 1.3 * pmax))
	if showlegend:
		axs[1, 0].legend([e10, i10, t10], ['e', 'i', 'tot'], fontsize=9)
	
	dvdt = dpdt / Ab / mp
	dt = 1.e-6  # timestep
	vtherm = np.sqrt(2. * ti * qe / Ab / mp)  # thermal ion velocity [m/s]
	t = [0]
	v_t = [vb]
	a_t = [np.interp(vb, vb_arr, dvdt)]
	while v_t[-1] > vtherm:
		t.append(t[-1] + dt)
		v_t.append(v_t[-1] - a_t[-1] * dt)
		a_t.append(np.interp(v_t[-1], vb_arr, dvdt))
	print('t_therm (numerical): {}s'.format(t[-1]))
	
	axs[0, 1].plot(t, np.array(v_t) / 1.e6, ls=ls)  # , label='Ti={}eV'.format(ti))
	axs[0, 1].set_ylabel('v(t) [Mm/s]')
	# axs[0, 1].legend()
	axs[1, 1].plot(t, np.array(a_t) / 1.e9, ls=ls)  # , label='Ti={}eV'.format(ti))
	axs[1, 1].set_ylabel('a(t) [Gm/s^2]')
	# axs[1, 1].legend()
	a = 1
	return axs


if __name__ == '__main__':
	
	demo = 1
	frac_demo = 0
	show1036170301 = 0
	if demo:
		te = 100  # plasma electron temp [eV]
		ti = 100  # plasma ion temp [eV]
		eb = 20000  # beam energy [eV]
		ni = ne = 1.e19  # plasma density [#/m^3]
		lnlamb = 15.  # unverified
		Zb = 1.  # beam charge number
		Ab = 1.  # beam atomic number  (?)
		axs = plot_collision_freqs()
		te = 200
		ti = 200
		plot_collision_freqs(axs=axs)
	if frac_demo:
		plot_fractional_energy_transfer(eb=11.e3, mb=1, mi=1, interps=[150, 250])
	if show1036170301:
		# run
		# psigs = ['\\te', '\\ti', '\\ni', '\\ne']
		# shot = 1036170301
		# for sig in psigs:
		# 	plot_sig(sig, shot)
		# in python directory to get the data reported below. values are eyeballed from the plots
		
		# EARLY
		te = 250.  # plasma electron temp [eV]
		ti = 50.  # plasma ion temp [eV]
		eb = 11000.  # beam energy [eV]
		ni = .75e19
		ne = 1.e19  # plasma density [#/m^3]
		lnlamb = 15.  # unverified
		Zb = 1.  # beam charge number
		Ab = 1.  # beam atomic number  (?)
		axs = plot_collision_freqs()
		# LATE
		ti = 250.
		ni = 1.75e19
		ne = 2.e19
		plot_collision_freqs(axs, ls='--', showlegend=False)
	
	plt.tight_layout()
	plt.show()
