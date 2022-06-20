import matplotlib.pyplot as plt
import numpy as np

# visualize beam footprint in machine

rtan = .21
phioffaxis = [0., 5 * np.pi / 180., 5 * np.pi / 180.]
phirotate = [0., np.pi / 2., 3. * np.pi / 2]
roffset = [0., .03, .03]

src_to_tangency = 2.57
focal_length = 1.8  # from NUBEAM setup
waist_to_tangency = src_to_tangency - focal_length
xwaist = 0
ywaist = np.sqrt(rtan ** 2 + waist_to_tangency ** 2)
zwaist = 0  # beam waist position [m] WLOG x=0
phibeamcenterline = np.arcsin(rtan / ywaist)
l = np.linspace(0, 2 * ywaist, int(np.ceil(200 * ywaist)))  # 1 pt per cm, transpose to keep correct format

rout, rin, theta = .67, .132, np.linspace(0, 2 * np.pi)
plt.plot(rout * np.cos(theta), rout * np.sin(theta), 'k')
plt.plot(rin * np.cos(theta), rin * np.sin(theta), 'k')

for i in [0, 1, 2]:
	roffsetxy = roffset[i] * np.sin(phirotate[i])
	z0 = zwaist + roffset[i] * np.cos(phirotate[i])
	x0 = xwaist - roffsetxy * np.cos(phibeamcenterline)  # negative from pos rotation giving offset towards HFS
	y0 = ywaist - roffsetxy * np.sin(phibeamcenterline)
	zpath = z0 + l * np.sin(phioffaxis[i]) * np.cos(phirotate[i])
	xpath = x0 + l * np.sin(phibeamcenterline - phioffaxis[i] * np.sin(phirotate[i]))
	ypath = y0 - l * np.cos(phibeamcenterline - phioffaxis[i] * np.sin(phirotate[i]))
	plt.plot(xpath, ypath)
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.tight_layout()
plt.axis('equal')
plt.show()

if __name__ == '__main__':
	pass
