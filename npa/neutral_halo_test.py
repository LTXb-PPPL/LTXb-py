import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

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
n0hh = neut_vars[
	'BOXN0H0'].data  # HIGHER GEN HALO NEUT DENS for all BOXES, N/CM**3 in format (NBBOX,NXBOX,NYBOX,NLBOX)'
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

beamtan = 21.3  # [cm] (tangency radius of beam: RTCENA in TR.DAT)
beamsrc_to_tan = 257  # [cm] (dist from beam source to tangency radius: XLBTNA in TR.DAT)
xsrc, ysrc = 0, np.sqrt(beamtan ** 2 + beamsrc_to_tan ** 2)
phi_src = np.arctan2(beamtan, beamsrc_to_tan)  # angle between src-machinecenter and beam-centerline

fig3, ax3 = plt.subplots(subplot_kw=dict(projection='polar'))
neut_halo.total_neut = neut_halo.boxn0[:,:,:,0] + neut_halo.boxn0h0[:,:,:,0] + neut_halo.boxn0hh[:,:,:,0]
n0tot = neut_halo.total_neut[:, :, :]
n0tot_mp = np.sum(n0tot, axis=1)  # sum over vertical coord
xbox, lbox = neut_halo.xbox, neut_halo.lbox
x_bb, y_bb = np.meshgrid(xbox, lbox)
r_bb, phi_bb = np.zeros_like(x_bb), np.zeros_like(x_bb)
for ix in np.arange(len(xbox)):
	for il in np.arange(len(lbox)):
		x_bb[il, ix] = xsrc + lbox[il] * np.sin(phi_src) - xbox[ix] * np.cos(phi_src)
		y_bb[il, ix] = ysrc - lbox[il] * np.cos(phi_src) - xbox[ix] * np.sin(phi_src)
		r_bb[il, ix] = np.sqrt(x_bb[il, ix] ** 2 + y_bb[il, ix] ** 2)
		phi_bb[il, ix] = np.arctan2(-x_bb[il, ix], y_bb[il, ix])  # neg here to get NB injecting co-Ip
# all neut stuff in [cm], convert to [m] to be consistent w/other plot
r_bb *= 1.e-2
cb3 = ax3.contourf(phi_bb, r_bb, n0tot_mp, levels=25)
cbar3 = fig3.colorbar(cb3)
cbar3.set_label('#/cm^3 (summed over vertical coord)')
ax3.set_rlim((0, .75))
ax3.set_theta_zero_location('N')
ax3.set_title('beam neutral halo')
a = 1

plt.tight_layout()
plt.show()

if __name__ == '__main__':
	pass
