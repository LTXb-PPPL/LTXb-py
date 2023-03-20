import numpy as np
import matplotlib.pyplot as plt

# controls
class WireControls:
	def __init__(self):
		self.nwire_bins = 100
		self.ngauss_bins = 100
		self.do_convergence_check = False
		
		# wire layout
		# wires placed at x-locations, span y-gap
		self.xgap, self.ygap = 24, 24  # [cm] aperture width in x, y (NOTE: 32 cm is TAE design)
		self.half_gap = self.xgap / 2.  # [cm] the aperture half-width (NOTE: TAE is 32cm width)
		self.nxwires, self.nywires = 7, 7
		
		# beam characteristics
		self.fwhmx = 8.  # cm
		self.fwhmy = 9.  # cm
		self.sigx = self.fwhmx / (2. * np.sqrt(2. * np.log(2)))  # [cm]
		self.sigy = self.fwhmy / (2. * np.sqrt(2. * np.log(2)))  # [cm]
		self.pbeam = 400.e3  # [W] total beam power
		self.tbeam = 10.e-3  # [s] beam pulse duration
		
		# wire characteristics: TAE design, 0.5mm diam wires 32 cm long
		self.w_diam = .05  # [cm]
		
		# 2d gaussian form- no covariance
		# z = A * np.exp(-(x - xw) ** 2 / (2. * sigx ** 2) - (y - yw) ** 2 / (2. * sigy ** 2))
		# so for given values of x, xw, sigx then z becomes 1d gaussian in y
		# total integral: V=2*pi*A*sigx*sigy=total beam power
		self.A = self.pbeam / (2. * np.pi * self.sigx * self.sigy)
		
		self.num_arrangements = 13

		# Tungsten characteristics
		self.c_mol = 24.27  # [J/mol/K] molar heat capacity
		self.rho = 19.25  # [g/cm^3]
		self.amu = 183.84  # [g/mol]
		self.c_sp = self.c_mol / self.amu  # [J/g/K] specific heat capacity
		self.k_con = 175.  # [W/m/K]  # thermal conductivity at 20C
		self.alpha = 4.3e-6  # [K^-1] thermal expansion coeff at 20C
