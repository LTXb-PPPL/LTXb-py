import matplotlib.pyplot as plt
import numpy as np
from op_scopes import nbi_ops

# Takeaways:
# 0 psi is optimal on GVC (better startup of Ibeam)
# -16 psi or so is optimal (flatter Ibeam)

gvc_psi = 511200 + np.array([162, 168, 80, 165, 126, 139, 153])  # psi: +4, +2, 0, -2, -4, -8, -12
shts = 511300 + np.array([65, 67])
nbi_ops(gvc_psi, arc_iv=False)
plt.show()
gva_psi = 511200 + np.array([52, 64, 75, 81, 88, 95])  # psi: 0, -4, -8, -12, -16, -20
nbi_ops(gva_psi, arc_iv=False)
plt.show()

clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']

if __name__ == '__main__':
	pass