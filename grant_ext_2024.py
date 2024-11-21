import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from bills_LTX_MDSplus_toolbox import *
import MDSplus

matplotlib.use('TkAgg')  # allows plotting in debug mode
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']

# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, sharex='col')
fig, (ax3, ax4) = plt.subplots(nrows=2, sharex='col')

def pull_transp_sig(conn, tshot, nodepath):
	conn.openTree(tree, tshot)
	node = conn.get(nodepath)
	data = node.data()
	dataunits = ' '.join(conn.get('units_of({})'.format(nodepath)).split())
	dim1 = conn.get('dim_of({}, 0)'.format(nodepath))
	dim1units = ' '.join(conn.get('units_of(dim_of({}, 0))'.format(nodepath)).split())
	return dim1, data


# nbi_ops([512528, 512534, 512560, 512565])  # 2 shots from 6/7/23 and 6/27/23 comparing 5 and 10 ms pulses
tree_shots = []  # [512528, 512534, 512560, 512565]
tree = 'ltx_nbi'
prefix = ''
for sh in tree_shots:
	t = get_tree_conn(sh, treename=tree)
	print(f'gathering data for shot {sh} occurred on {get_data(t, ".metadata:timestamp")}')
	(tibeam, ibeam) = get_data(t, f'{prefix}.source_diags.i_hvps')
	(tbeam, vbeam) = get_data(t, f'{prefix}.source_diags.v_hvps')
	ibeam = np.interp(tbeam, tibeam, ibeam)
	pbeam = ibeam * vbeam
	ax1.plot(tbeam, pbeam, label=sh)
	ax2.plot(tbeam, vbeam)

# transp
tshots = [1088901601]#, 1088902202]
tree = 'transp_ltx'
conn_name = 'transpgrid'
conn_name2 = 'transpgrid2.pppl.gov'
conn = MDSplus.Connection(conn_name)
for tshot in tshots:
	print(f'gathering transp data for shot {tshot}')
	ptime, pbeam = pull_transp_sig(conn, tshot, '\\pinj')
	vtime, vbeam = pull_transp_sig(conn, tshot, '\\einj')
	
	ax3.plot(ptime, pbeam, label=tshot)
	ax4.plot(vtime, vbeam)

# ax1.set_ylabel('pbeam')
# ax2.set_ylabel('vbeam')
ax3.set_ylabel('transp pbeam')
ax4.set_ylabel('transp vbeam')
ax3.legend()
plt.show()

if __name__ == '__main__':
	pass
