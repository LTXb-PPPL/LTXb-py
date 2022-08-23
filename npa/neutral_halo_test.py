import matplotlib.pyplot as plt
import numpy as np
from transp_code.transp_classes import Halo3D

n0_fn = '//samba/wcapecch/transp/t103617/103617C01_boxn0_4.cdf'
neut_halo = Halo3D(n0_fn)
neut_vars = neut_halo.vars

xb = neut_vars['XBOX'].data  # X grid in 3D BOX, CM
yb = neut_vars['YBOX'].data  # Y grid in 3D BOX, CM
lb = neut_vars['LBOX'].data  # L grid in 3D BOX, CM

# dimensions below: ['NLBOX', 'NYBOX', 'NXBOX', 'NBBOX']
n0 = neut_vars['BOXN0'].data  # BEAM NEUTRAL DENSITIES for all BOXES, N/CM**3 in format (NBBOX,NXBOX,NYBOX,NLBOX)
n0h0 = neut_vars['BOXN0H0'].data  # 0 GEN HALO NEUTRAL DENS for all BOXES, N/CM**3 in format (NBBOX,NXBOX,NYBOX,NLBOX)'
n0hh = neut_vars['BOXN0H0'].data  # HIGHER GEN HALO NEUT DENS for all BOXES, N/CM**3 in format (NBBOX,NXBOX,NYBOX,NLBOX)'
n0tot = n0 + n0h0 + n0hh  # beam neutrals plus 0 gen plus higher gen neutrals

fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, figsize=(18, 10))
cont1 = ax1.contourf(xb, lb, np.sum(n0tot[:, :, :, 0], axis=1), levels=25)
cont2 = ax2.contourf(xb, lb, np.sum(n0[:, :, :, 0], axis=1), levels=cont1.levels)
cont3 = ax3.contourf(xb, lb, np.sum(n0h0[:, :, :, 0], axis=1), levels=cont1.levels)
cont4 = ax4.contourf(xb, lb, np.sum(n0hh[:, :, :, 0], axis=1), levels=cont1.levels)
ax1.set_title('total neutral density')
ax2.set_title('beam neutral density')
ax3.set_title('0 gen neutrals')
ax4.set_title('higher gen neutrals')

cbar = fig.colorbar(cont1)
cbar.ax.set_ylabel('sort of N/CM^3')
ax1.set_ylabel('L (cm)')
ax2.set_xlabel('X (cm)')

plt.tight_layout()
plt.show()


if __name__ == '__main__':
	pass
