import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from bills_LTX_MDSplus_toolbox import *
import MDSplus

matplotlib.use('TkAgg')  # allows plotting in debug mode
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']


class PullSignal:
	"""A signal class to provide an organizational scheme for storing
	timeseries data retrieved from an MDSplus signal node. To collect
	the data from an MDSplus signal node, provide the shot number,
	nodepath, and tree.
	"""
	
	def __init__(self, shot, nodepath, tree='transp_ltx'):
		"""Instatiating a signal will make an attempt to collect data
		for a particular node on the MDSplus server on
		dave.physics.wisc.edu.  Specify an integer for the shot number
		(e.g. 1161127010) and a string that describes the node path
		(e.g. '\ip'). The default tree is 'mst', but you can choose
		others available on dave.physics.wisc.edu."""
		if tree == 'mst':
			conn_name = 'dave.physics.wisc.edu'
			conn_name2 = None
			goodshot_node_test = '\ip'
		elif tree == 'ltx_b':
			conn_name = 'lithos'
			conn_name2 = None
			goodshot_node_test = '.diagnostics.magnetic.ip_rog'
		elif tree == 'transp_ltx':
			conn_name = 'transpgrid'
			conn_name2 = 'transpgrid2.pppl.gov'
			goodshot_node_test = '\\dvol'
		
		# Initialize the fields of the simpleSignal class to be the input of
		# the constructor.
		self.shot = shot
		self.nodepath = nodepath
		self.tree = tree
		self.data = None
		self.dataunits = None
		self.ndim = None
		self.dim1 = None
		self.dim1units = None
		self.dim2 = None
		self.dim2units = None
		self.dim3 = None
		self.dim3units = None
		self.error = None
		self.name = None
		self.label = None
		self.goodshot = None  # indicates whether or not a shot exists in tree with any data
		
		# Get name of signal, last name of the nodepath
		self.name = self.nodepath.split(':')[-1].split('.')[-1]
		
		# Look at standard signal as indicator of whether or not data for this shot exists at all
		try:
			conn = MDSplus.Connection(conn_name)
			conn.openTree(tree, shot)
			testfordata = conn.get(goodshot_node_test)
			self.goodshot = True
			# Tree could have more data associated with a Signal, attempt
			# to get this information, but deal with the possibility that
			# this data may not exist for a given signal.
			try:
				node = conn.get(nodepath)
				self.data = node.data()
				self.dataunits = ' '.join(conn.get('units_of({})'.format(nodepath)).split())
				self.ndim = self.data.ndim
				self.dim1 = conn.get('dim_of({}, 0)'.format(nodepath))
				self.dim1units = ' '.join(conn.get('units_of(dim_of({}, 0))'.format(nodepath)).split())
				if self.ndim > 1:  # 2 or 3 dimensional data, get data/units for 2nd dimension
					self.dim2 = conn.get('dim_of({}, 1)'.format(nodepath))
					self.dim2units = ' '.join(conn.get('units_of(dim_of({}, 1))'.format(nodepath)).split())
				if self.ndim > 2:  # 3 dimensional data
					self.dim3 = conn.get('dim_of({}, 2)'.format(nodepath))
					self.dim3units = ' '.join(conn.get('units_of(dim_of({}, 2))'.format(nodepath)).split())
			except:
				print("{0} is not available for shot {1}".format(self.name, self.shot))
				return
		except:
			if conn_name2 is not None:
				try:
					conn = MDSplus.Connection(conn_name2)
					conn.openTree(tree, shot)
					testfordata = conn.get(goodshot_node_test)
					self.goodshot = True
					# Tree could have more data associated with a Signal, attempt
					# to get this information, but deal with the possibility that
					# this data may not exist for a given signal.
					try:
						node = conn.get(nodepath)
						self.data = node.data()
						self.dataunits = ' '.join(conn.get('units_of({})'.format(nodepath)).split())
						self.ndim = self.data.ndim
						self.dim1 = conn.get('dim_of({}, 0)'.format(nodepath))
						self.dim1units = ' '.join(conn.get('units_of(dim_of({}, 0))'.format(nodepath)).split())
						if self.ndim > 1:  # 2 or 3 dimensional data, get data/units for 2nd dimension
							self.dim2 = conn.get('dim_of({}, 1)'.format(nodepath))
							self.dim2units = ' '.join(conn.get('units_of(dim_of({}, 1))'.format(nodepath)).split())
						if self.ndim > 2:  # 3 dimensional data
							self.dim3 = conn.get('dim_of({}, 2)'.format(nodepath))
							self.dim3units = ' '.join(conn.get('units_of(dim_of({}, 2))'.format(nodepath)).split())
					except:
						print("{0} is not available for shot {1}".format(self.name, self.shot))
						return
				except:
					self.goodshot = False
					return
			else:
				self.goodshot = False
				return
		
		# Let's see if there's any more meta data like units and error bars.
		
		try:
			self.label = ' '.join(conn.get('{}.label'.format(nodepath)).split())
		except:
			pass
		if self.label is None:
			self.label = self.name
		
		# Close connection to prevent hanging processes.
		conn.closeAllTrees()
	
	def plot(self, title=None, contour=False, label=None, ax=None):
		"""Plot the signal vs. time using as much of the stored
		 information as is available to annotate the figure. Set the
		 color of the plot line using, e.g. color='b' for blue."""
		
		if ax is None:
			fig = plt.figure()
		if self.ndim == 1:
			plt.plot(self.dim1, self.data, label=label)
			plt.xlabel('x ({})'.format(self.dim1units))
			plt.ylabel('{} ({})'.format(self.label, self.dataunits))
		elif self.ndim == 2:
			if contour:
				if ax is None:
					ax = plt.axes()
				ax.contourf(self.dim1, self.dim2, self.data)
				plt.title('{} ({})'.format(self.label, self.dataunits))
			else:
				if ax is None:
					ax = plt.axes(projection='3d')
				xx, yy = np.meshgrid(self.dim1, self.dim2)
				ax.plot_surface(xx, yy, self.data, cmap='viridis')
				ax.set_zlabel('{} ({})'.format(self.label, self.dataunits))
			ax.set_xlabel('x ({})'.format(self.dim1units))
			ax.set_ylabel('y ({})'.format(self.dim2units))
			ax.set_title(title)
		elif self.ndim == 3:
			num = len(self.data)  # len here gives length of 1st dimension
			xx, yy = np.meshgrid(self.dim2, self.dim3)  # assume dim1 is time axis
			iplot = [int(ii) for ii in np.linspace(0, len(self.dim1) - 1, num=3)]
			for n, i in enumerate(iplot):
				if contour:
					ax = fig.add_subplot(3, 1, n + 1)
					ax.contourf(self.dim2, self.dim3, self.data[:, :, i],
					            label='{} ({})'.format(self.label, self.dataunits))
					ax.legend()
				else:
					ax = fig.add_subplot(3, 1, n + 1, projection='3d')
					ax.plot_surface(xx, yy, self.data[:, :, i], cmap='viridis')
					ax.set_zlabel('{} ({})'.format(self.label, self.dataunits))
					ax.set_title('t: {}'.format(self.dim1[i]))
				ax.set_xlabel('x ({})'.format(self.dim1units))
				ax.set_ylabel('y ({})'.format(self.dim2units))


# transp
def plot_timevar_comparison():
	# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, sharex='col')
	# ax2r, ax4r = ax2.twinx(), ax4.twinx()
	fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True)
	ax2r = ax2.twinx()
	tshots = [1088901601, 1078011601, 1088901602, 1078011602]
	lbls = ['10 ms exp', '5 ms exp', '10 ms avg power', '5 ms avg power']
	tree = 'transp_ltx'
	conn_name = 'transpgrid'
	conn_name2 = 'transpgrid2.pppl.gov'
	conn = MDSplus.Connection(conn_name)
	for i, tshot in enumerate(tshots):
		print(f'gathering transp data for shot {tshot}')
		pbeam = PullSignal(tshot, '\\pinj')  # [W]
		bpte = PullSignal(tshot, '\\bpte')  # [W]
		
		ax1.plot(pbeam.dim1 * 1000., pbeam.data / 1000., label=lbls[i])
		ax2r.plot(bpte.dim1 * 1000., bpte.data / 1000.)
	
	ax1.set_ylabel('pbeam (kW)')
	ax2r.set_ylabel('elec heating (kW)')
	for ax in (ax1, ax2):
		ax.set_xlabel('time (ms)')
		ax.set_xlim((460, 485))
	ax1.legend()
	plt.tight_layout()
	plt.show()


if __name__ == '__main__':
	plot_timevar_comparison()
