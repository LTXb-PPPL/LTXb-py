import matplotlib.pyplot as plt
import numpy as np
from nptdms import TdmsFile

dir = 'C:/Users/wcapecch/Dropbox/work_stuff/nbi/example_shot_data/'
for fn in [f'{dir}adc75index.tdms']:
	tdmsfile = TdmsFile.read(fn)
	df = tdmsfile.as_dataframe()
	# for i in np.arange(30):
	# 	plt.plot(df[df.columns[2 * i]], df[df.columns[2 * i + 1]], label=f'{2*i}')
	# plt.legend()
	# plt.show()
	a = 1

if __name__ == '__main__':
	pass
