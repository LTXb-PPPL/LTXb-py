import pickle

import matplotlib.pyplot as plt
import numpy as np
from toolbox.particle_trajectory_tools import get_generic_path
from toolbox.helpful_stuff import ltx_above, ltx_limiter

clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']

nth = 3
fsav2d = f'Z:/PycharmProjects/LTXb-py/neutral_beam/neut_halo_fueling_ds{nth}.pkl'
a = np.linspace(1, 5)
pickle.dump(a, open(fsav2d, 'wb'))
print(f'saved: {fsav2d}')
a = pickle.load(open(fsav2d, 'rb'))
plt.plot(a)
plt.show()

if __name__ == '__main__':
	pass
