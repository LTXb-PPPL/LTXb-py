import matplotlib.pyplot as plt
import numpy as np
from bills_LTX_MDSplus_toolbox import *

# check if power in arc signals roughly corresponds to drop in voltage on arcPS

shot = 107977  # 10ms pulse @ 13kV
if shot > 200000:
	tree = 'ltx_nbi'
	nbi_only = True
	prefix = ''
else:
	tree = 'ltx_b'
	nbi_only = False
	prefix = '.oper_diags.ltx_nbi'
t = get_tree_conn(shot, treename=tree)
(tiarc, iarc) = get_data(t, f'{prefix}.source_diags.i_arc')
(tarc, varc) = get_data(t, f'{prefix}.source_diags.v_arc')
iarc = np.interp(tarc, tiarc, iarc)

# correct offset
tstrt = np.where(tarc < .464)
iarc -= np.mean(iarc[tstrt])
varc -= np.mean(varc[tstrt])

parc = iarc * varc
arc_energy_J = np.sum(parc[1:] * (tarc[1:] - tarc[:-1]))

# arc PS bank
c = 1.e4 * 1.e-6  # 10,000 uF caps
num_caps = 16  # 16 caps in parallel
vset = 330  # ps saturates to around 330 V
ps_stored_energy = 0.5 * c * num_caps * vset ** 2
ufinal = ps_stored_energy - arc_energy_J
vfinal = np.sqrt(ufinal * 2 / (c * num_caps))
# pred_vdrop = np.sqrt(arc_energy_J * 2 / (c * num_caps))

print(f'energy stored: {ps_stored_energy:.2f} J')
print(f'energy in arc: {arc_energy_J:.2f} J')
print(f'predicted final voltage: {vfinal:.2f} V')
# print(f'voltage predicted to drop to: {vset - vfinal:.2f} V')

clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
fig, ax = plt.subplots()
axr = ax.twinx()
ax.plot(tarc, iarc, label='iarc')
axr.plot(tarc, varc, label='varc', c=clrs[1])
ax.set_ylabel('I_arc', c=clrs[0])
axr.set_ylabel('V_arc', color=clrs[1])
ax.set_xlabel('t (s)')
ax.set_ylim(bottom=0)
axr.set_ylim(bottom=0)

plt.show()

if __name__ == '__main__':
	pass
