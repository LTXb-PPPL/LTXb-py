from toolbox.helpful_stuff import SimpleSignal
from toolbox.ionization_cross_sections import *

"""
Idea here- at our injection energy, only 25-30% of beam is fueling (via ion- or electron-impact ionization).
The CX fraction dominates and simply swaps a thermal ion & fast neutral to a thermal neutral & fast ion.
The 70-75% of the beam that creates a population of thermal neutrals still has a chance to fuel via secondary (and
tertiary...) ionization processes and could contribute to beam fueling beyond the 25-30% of the initial ionization processes.
"""

# shot = 1074580301
shot = 1088900311  # 0302, 0310, 0311
pvol = SimpleSignal(shot, '\\pvol')
pvolf = SimpleSignal(shot, '\\pvolf')
einj = SimpleSignal(shot, '\\einj')  # beam injection energy (eV)
dvol = SimpleSignal(shot, '\\dvol')  # zone volume (cm^-3)
sbtot = SimpleSignal(shot, '\\sbtot')  # total ion source (beam+halo) (num/sec/cm^-3)  VS x"r/a" and time
sinj = SimpleSignal(shot, '\\sinj')  # beam neutrals injected (num/sec)
shin = SimpleSignal(shot, '\\sbshine_h')  # shine-through (num/sec)
ionsrc = np.sum(sbtot.data * dvol.data, axis=1)

dt = sinj.dim1[1:] - sinj.dim1[:-1]
ninj = np.sum(sinj.data[1:] * dt)
nshine = np.sum(shin.data[1:] * dt)
nsrc = np.sum(ionsrc[1:] * dt)

te = 150.  # estimate of elec temp
eb = max(
	einj.data)  # NOTE that for run 107458C01 at least, EINJ is a very weird signal. Nonzero before beam, doesn't end at beam shutoff
# not_shinethrough =
sigv_ii = tabulated_iimpact_ionization(eb, 1.)
sigv_ie = tabulated_eimpact_ionization(eb, te)
sigv_cx = tabulated_i_ch_ex_ionization(eb, 1.)
fuel_frac_est = (sigv_ii + sigv_ie) / (sigv_cx + sigv_ii + sigv_ie) * (1. - nshine / ninj)
print(f'shinethrough: {nshine / ninj * 100:.2f}%')

fig, ax = plt.subplots()
p1, = ax.plot(sinj.dim1, sinj.data, label=f'injec: {ninj:.2e}')
p2, = ax.plot(sinj.dim1, ionsrc, label=f'ion src: {nsrc:.2e}')
p3, = ax.plot(np.nan, np.nan, label=f'fueling: {nsrc / ninj * 100:.1f}%')
p4, = ax.plot(np.nan, np.nan, label=f'1st pass fueling est: {fuel_frac_est * 100:.1f}%\n (assuming Te={te} eV')
leg1 = ax.legend(handles=[p1, p2], loc='upper left', title='#/sec')
ax.legend(handles=[p3, p4], loc='center left')
ax.add_artist(leg1)
plt.title(f'{shot}')

plt.show()

# l1, = ax1.plot(np.nan, np.nan, 'k-', label='$I_p$ (kA)')
# l2, = ax1.plot(np.nan, np.nan, 'k--', label='$P_{OH} (kW)$')
# leg1 = ax1.legend(handles=[l1, l2], loc='upper right', labelcolor='linecolor', frameon=False, fontsize=fs2)
# l3, = ax1.plot(np.nan, np.nan, 'k', label='No Beam', ls='none')
# l4, = ax1.plot(np.nan, np.nan, 'r', label='With Beam', ls='none')
# ax1.legend(handles=[l3, l4], loc='upper left', labelcolor='linecolor', frameon=False, fontsize=fs2, handlelength=0)
# ax1.add_artist(leg1)
a = 1
