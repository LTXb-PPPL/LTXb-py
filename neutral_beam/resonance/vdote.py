import matplotlib.pyplot as plt
import numpy as np
from helpful_stuff import read_eqdsk
import pickle
from scipy import integrate
from scipy.interpolate import interp2d


def vdote(mnum, nnum, fepm, psimode, phi0, verbose=False, plot=False, recon=None, orbit=None):
	singular = False
	if len(fepm) == 1 and len(psimode) == 1 and len(phi0) == 1:
		singular = True
	
	nf, nr, n0 = 20., 20., 1.
	# out = read_eqdsk(recon)
	orbit_dat = pickle.load(open(orbit, 'rb'))
	orb = orbit_dat['orb']
	out = orb['OUT']
	t = np.transpose(orb['T_ORB'])  # default orbit has 20,000 pts from 0-66.4microsec
	
	# TODO; do we need all these transposes below??
	# COMPUTE PARTICLE vr ARRAY (notebook2 pg 135)
	r_mag = orb['OUT']['R0']  # [m]
	Rmaj_par = np.sqrt(orb['X'] ** 2 + orb['Y'] ** 2)
	psifunc = interp2d(out['X'], out['Y'], out['NORM_PSIXY'],
	                   kind='cubic')  # function to return interpolated psi values
	normpsi_par = np.copy(Rmaj_par) * 0
	for ipsipar in np.arange(len(normpsi_par)):
		normpsi_par[ipsipar] = psifunc(Rmaj_par[ipsipar], orb['Z'][ipsipar])[
			0]  # have to do pt by pt otherwise it grids R, z inputs
	rminor_par = np.sqrt((Rmaj_par - r_mag) ** 2 + orb['Z'] ** 2)
	par_theta = np.arctan2(orb['Z'], Rmaj_par - r_mag)
	par_theta = np.transpose((par_theta + 2 * np.pi) % (2 * np.pi))  # on interval [0, 2pi)
	par_phi = np.transpose((orb['PHI'] + 2 * np.pi) % (2 * np.pi))
	par_vr = np.transpose(orb['VR'] * np.cos(par_theta) + orb['VZ'] * np.sin(par_theta))
	par_vtheta = np.transpose(orb['VR'] * np.sin(par_theta) - orb['VZ'] * np.cos(par_theta))
	par_vphi = np.transpose(orb['VPHI'])
	
	if psimode is None:
		psimode = np.mean(normpsi_par)
	
	# here we solve for the exact resonant frequency based on orbit
	orbit_freq = abs(nnum * orb['CIQ']['W_PHI'] - mnum * orb['CIQ']['W_THETA']) / 2. / np.pi
	if fepm is None:
		fepm = [orbit_freq]
	
	power = t * 0.
	vdote_arr = np.zeros((len(psimode), len(fepm), len(phi0)))
	trend = np.copy(vdote_arr)
	
	last = 0
	npts = len(psimode) * len(fepm) * len(phi0)
	i_mp = np.where(abs(out['Y']) == min(abs(out['Y'])))[0][0]  # index closest to midplane (should be = 32)
	for rr in range(len(psimode)):
		if verbose:
			print(len(psimode) - rr)
		for ff in range(len(fepm)):
			for oo in range(len(phi0)):
				perc = (oo + ff * len(phi0) + len(fepm) * len(phi0) * rr) * 100. / (
						len(psimode) * len(fepm) * len(phi0))
				if perc - last > 10:
					print(f'{np.floor(perc)}% orbit done')
					last = np.floor(perc)
				
				psi_v = psimode[rr]
				fepm_v = fepm[ff]
				
				# TODO: This might not work- field line pitch isn't constant on flux surface, mode structure is complicated
				# todo: defining mode rotation based on outboard midplane pitch- how wrong is that?
				x_outboardmidplane = max(out['X'][np.where(out['NORM_PSIXY'][i_mp, :] <= psimode[rr])])
				b_theta = np.interp(x_outboardmidplane, out['X'], out['BP'])
				b_phi = np.interp(x_outboardmidplane, out['X'], out['BPHI'])
				
				# look at variation of q on flux surface
				if 0:
					cs = plt.contour(out['X_XY'], out['Y_XY'], out['NORM_PSIXY'], levels=[psimode[rr]])
					for item in cs.collections:
						for i in item.get_paths():
							v = i.vertices  # vertices of contour
					xcont, ycont = v[:, 0], v[:, 1]
					bp_func = interp2d(out['X'], out['Y'], out['BP_XY'], kind='cubic')
					bphi_func = interp2d(out['X'], out['Y'], out['BPHI_XY'], kind='cubic')
					b_theta = bp_func(max(xcont), 0.)  # value of bp on outboard midplane
					b_phi = bphi_func(max(xcont), 0.)  # value of bphi on outboard midplane
					qfunc = interp2d(out['X'], out['Y'], out['Q_XY'], kind='cubic')
					q_fluxsurf, bp_fs, bphi_fs = np.copy(xcont) * 0, np.copy(xcont) * 0, np.copy(xcont) * 0
					for ifs in np.arange(len(q_fluxsurf)):
						q_fluxsurf[ifs] = qfunc(xcont[ifs], ycont[ifs])[0]
						bp_fs[ifs] = bp_func(xcont[ifs], ycont[ifs])[0]
						bphi_fs[ifs] = bphi_func(xcont[ifs], ycont[ifs])[0]
				
				# original code on next 3 lines
				# b_theta = np.interp(out['bp'], out['rhov'], psimode[rr])
				# b_phi = np.interp(out['bphi'], out['rhov'], psimode[rr])
				# q = np.interp(out['q'], out['rhov'], psimode[rr])
				
				# phi = np.linspace(0, 2 * np.pi, endpoint=True)  # findgen(1001)/1000.*2*!pi
				# theta = np.copy(phi)  # findgen(1001)/1000.*2*!pi
				
				amp0 = 1.  # arbitrary anyway
				sigma = .02  # guassian radial spread of mode
				
				# IF 0 THEN BEGIN ;as an aside, plot br contours
				# 	 rp = (findgen(1001)/1000.-.5)*1.04+1.5
				# 	 zp = (findgen(1001)/1000.-.5)*1.04
				#
				# 	 epm_pol = make_array(n_elements(rp),n_elements(zp))
				# 	 FOR r=0,n_elements(rp)-1 DO BEGIN
				# 	 FOR z=0,n_elements(zp)-1 DO $
				# 		epm_pol[r,z] = amp0*exp(-(sqrt((rp[r]-1.5)^2+zp[z]^2) - psimode) ^ 2 $
				# 					   /2./sigma^2)*sin(m*atan(zp[z],rp[r]-1.5)+!pi/3.)
				# 	ENDFOR
				# 	 th = findgen(1001)/1000.*!pi
				# 	 epm_tor = make_array(n_elements(th),n_elements(phi))
				# 	 FOR p=0,n_elements(phi)-1 DO BEGIN
				# 		  FOR t=0,n_elements(th)-1 DO BEGIN
				# 			   epm_tor[t,p] = amp0*sin(m*th[t])*sin(n*phi[p])
				# 			  ENDFOR
				# 		 ENDFOR
				#
				# 	 putwindow,5
				# 	 !p.multi=[0,2,1]
				# 	 nl=50
				# 	 lev=findgen(nl)/(nl-1)*(max(epm_pol)-min(epm_pol))+min(epm_pol)
				# 	 contour,epm_pol,rp,zp,xtit='radius (m)',ytit='z (m)', $
				# 			 lev=lev,/iso,/fill
				# 	 contour,epm_tor,th/!pi,phi/!pi,levels=findgen(25)/24.*(max(epm_tor)- $
				# 			 min(epm_tor))+min(epm_tor),/fill,xr=[0,1.],yr=[0,2.], $
				# 			 xtit='theta/pi',ytit='phi/pi'
				#
				# ENDIF ;end of plot br contours
				
				# field_line = (phi / q) % (2 * np.pi)
				#
				# Rphi = np.linspace(0, 2 * np.pi * r_mag)  # findgen(1001)/1000. * 2*!pi*r_maj
				# rtheta = np.linspace(0, 2 * np.pi * psi_v)  # findgen(1001)/1000. * 2*!pi*r_min_v
				
				omega = 2 * np.pi * fepm_v
				delta_phi = omega * t + phi0[oo]
				line_pitch = b_theta / b_phi
				delta_theta = delta_phi * line_pitch
				
				# have delta_theta,delta_phi of mode as functions of time
				# as well as par_theta,par_phi location of particle as psifunc of time
				
				amp = amp0 * np.exp(-(rminor_par - psimode[rr]) ** 2 / 2. / sigma ** 2)
				epm = amp * np.sin(nnum * par_phi + delta_phi) * np.sin(mnum * par_theta + delta_theta)
				
				# if 0:  # compute epm amp at particle location as psifunc of time
				# 	power = epm*par_vr
				# 	vdote[rr,ff,oo] = int_tabulated(t,epm*par_vr)
				# compute simplified vdote using curl (pg. 136)
				etheta = -1. / fepm[ff] * amp * np.cos(nnum * par_phi + delta_phi) * np.sin(
					mnum * par_theta + delta_theta)
				ephi = 1. / fepm[ff] * amp * np.sin(nnum * par_phi + delta_phi) * np.cos(
					mnum * par_theta + delta_theta)
				power = etheta * par_vtheta + ephi * par_vphi
				vdote_arr[rr, ff, oo] = integrate.simpson(power, t)
				int_power = np.cumsum(power * (t[1] - t[0]))  # assumes regular intervals in t array
				fit = np.polyfit(t, int_power, 1)
				# track slope of total power transfer (ORIGINAL USED abs(fit[0])- don't want that! want to see if power is going to particle or to mode)
				trend[rr, ff, oo] = fit[0]
			
			# end of offset array
	# end of fepm array
	# end of psimode array
	
	imax = np.where(trend.flatten() == max(trend.flatten()))[0][0]
	omax = int(imax / (len(psimode) * len(fepm)))
	imax2 = int(imax % (len(psimode) * len(fepm)))
	# row index (pmax) = int(imax2/len(row)) where len(row) = len(fepm)
	# col index (fmax) = imax2 % len(row)
	pmax = int(imax2 / len(fepm))
	fmax = int(imax2 % len(fepm))
	
	if not singular:
		optimum = vdote(mnum, nnum, [fepm[fmax]], [psimode[pmax]], [phi0[omax]], verbose=False, plot=False, recon=recon,
		                orbit=orbit)
		t = optimum['t']
		power = optimum['power']
		int_power = optimum['int_power']
	
	resonance = {'vdote': vdote_arr, 'psimode': psimode, 'fepm': fepm, 'phi0': phi0, 'fmax': fmax, 'pmax': pmax,
	             'omax': omax, 't': t, 'power': power, 'int_power': int_power, 'trend': trend}
	
	# PLOTTING
	# vdote_arr *= 1.e5  # arbitrary scale anyway
	if plot and len(fepm) > 1 and len(psimode) > 1:
		fig, (ax1, ax2) = plt.subplots(ncols=2)
		ax1.plot(out['RVES'], out['ZVES'])
		ax1.plot(orb['R'][:2000], orb['Z'][:2000])
		ax1.contour(out['X_XY'], out['Y_XY'], out['NORM_PSIXY'], levels=[psimode[pmax]], zorder=5)
		
		nl = 50
		levs = np.linspace(min(trend.flatten()), max(trend.flatten()), num=nl, endpoint=True)
		ax2.contourf(psimode, fepm / 1000., trend[:, :, omax].transpose(), levels=levs)
		
		ax1.set_xlabel('X')
		ax1.set_ylabel('Y')
		ax2.set_xlabel('psi_mode')
		ax2.set_ylabel('f (kHz)')
		ax2.set_title('vdote (arb)')
		ax2.plot(psimode[pmax], fepm[fmax] / 1000., 'rx')
		
		if verbose:
			print(
				f'resonant at:\nr_mode = {psimode[pmax] * 100.}cm\nfepm = {fepm[fmax] / 1000.}kHz\nphi_0 = {phi0[omax]}')
	# todo: sort out the rest below
	# !p.multi=[0,2,1]
	# yt = 6
	# ytn = replicate(' ',yt+1)
	# plotval = 1.1*[transpose(vdote[rmax,*,omax]),vdote[*,fmax,omax]]
	# yr = [min(plotval),max(plotval)]
	# plot,fepm/1000.,vdote[rmax,*,omax],psym=4,xtit='fepm (kHz)',yr=yr,/yst,$
	# 	 ytit='vdote (arb)',charsize=cs,yticks=yt,xmargin=[10,-3]
	# hline,0
	# plots,fepm[fmax]/1000.,vdote[rmax,fmax,omax],psym=4,color=!red,thick=2
	# billegend,['r_mode = ' + trim(psimode[rmax] * 100.)],/ right, box= 0, chars=cs
	# ;vdote vs r_mode
	# plot, psimode * 100., vdote[*, fmax, omax],psym= 4, xtit= 'r_mode (cm)', yr= yr,/ yst, $
	# 	 chars=cs,yticks=yt,ytickn=ytn,xmargin=[4,3]
	# hline,0
	# plots, psimode[rmax] * 100., vdote[rmax, fmax, omax], psym= 4, color=!red, thick=2
	# billegend,['fepm = '+trim(fepm[fmax]/1000.)],/right,box=0,chars=cs
	#
	# putwindow,2
	# !p.multi=[0,1,2]
	# xt = 8
	# xtn = replicate(' ',xt+1)
	# t = t*1.e6 ;convert to microseconds
	# int_t = int_t*1.e6
	# plot,t,power,ytit='power(t) (arb)', $
	# 	 chars=cs,xticks=xt,xtickn=xtn,ymar=[2,2]
	# int_power = int_power*1.e5 ;arbitrary anyway
	# plot,int_t,int_power,ymar=[4,-1],xticks=xt,chars=cs, $
	# 	 ytit=textoidl('int(P)'),xtit=textoidl('time (\mu s)')
	#
	# putwindow,4
	# !p.multi=[0,2,1]
	# plot,data.r[0:10000],data.z[0:10000],/iso,xtit='R (m)',ytit='Z (m)'
	# plot,data.x[0:10000],data.y[0:10000],/iso,xtit='X (m)',ytit='Y (m)'
	#
	# IF keyword_set(verbose) THEN BEGIN
	# 	 print,'calculated rmode = '+trim(mean(rminor_par)*100.)+'cm'
	# 	 print,'calculated orbit freq = '+trim(orbit_freq/1000.)+'kHz'
	# 	ENDIF
	# ENDIF ;end of plotting
	
	return resonance


if __name__ == '__main__':
	pass
