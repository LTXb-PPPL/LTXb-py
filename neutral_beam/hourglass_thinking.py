import matplotlib.pyplot as plt
import numpy as np

# IR imaging of the beam dump shows an hourglass-like emission
# I'm wondering if that could be a gaussian-like beam spot overlaid on a vertical gradient
# vertical gradient could be plasma erosion of the dump surface due to stray field in the gap (?)

x = np.linspace(-50, 50)
y = np.linspace(-50, 50)
xx, yy = np.meshgrid(x, y)

gsig = 100
lsig = 20
center_frac = .8


def gauss(x, y):
	return np.exp(-(x / gsig) ** 2 - (y / gsig) ** 2)


def liner(y):
	return 1 - (1 - center_frac) * np.exp(-(y / lsig) ** 2)


zgaus = gauss(xx, yy)
zlin = liner(yy)

fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(15, 5), tight_layout=True)
ax1.contourf(xx, yy, zlin)
ax2.contourf(xx, yy, zgaus)
ax3.contourf(xx, yy, zlin * zgaus)

fig2, (aa, bb, cc, dd, ee, ff) = plt.subplots(ncols=6, sharey=True, figsize=(14, 3), tight_layout=True)
center_frac = 1.
for ax in [aa, bb, cc, dd, ee, ff]:
	ax.contourf(xx, yy, zgaus * liner(yy))
	ax.set_title(f'{center_frac:.2f}')
	center_frac -= .2

for ax in [ax1, ax2, ax3, aa, bb, cc, dd, ee, ff]:
	ax.axis('off')
	ax.axis('equal')

ax1.set_title('dump emissivity', fontsize=16)
ax2.set_title('beam power', fontsize=16)
ax3.set_title('thermal beam emission', fontsize=16)
fig2.suptitle('thermal beam emission\nvs central normalized emissivity')

plt.show()
a = 1
