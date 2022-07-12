import matplotlib.pyplot as plt
import numpy as np
from helpful_stuff import read_eqdsk
import pickle
from scipy import integrate


def vdote(mnum, nnum, fepm, psimode, phi0, verbose=False, plot=False, recon=None, orbit=None, nosave=False,
          finescan=False):
	singular = False
	if len(fepm) == 1 and len(psimode) == 1 and len(phi0) == 1:
		singular = True
	
	nf, nr, n0 = 20., 20., 1.
	out = read_eqdsk(recon)
	orbit_dat = pickle.load(open(orbit, 'rb'))
	orb = orbit_dat['orb']
	t = np.transpose(orb['T_ORB'])  # default orbit has 20,000 pts from 0-66.4microsec
	
	# TODO; do we need all these transposes below??
	# COMPUTE PARTICLE vr ARRAY (notebook2 pg 135)
	r_mag = orb['OUT']['R0']  # [m]
	Rmaj_par = np.sqrt(orb['X'] ** 2 + orb['Y'] ** 2)
	rminor = np.sqrt((Rmaj_par - r_mag) ** 2 + orb['Z'] ** 2)
	par_theta = np.arctan2(orb['Z'], Rmaj_par - r_mag)
	par_theta = np.transpose((par_theta + 2 * np.pi) % (2 * np.pi))  # on interval [0, 2pi)
	par_vr = np.transpose(orb['VR'] * np.cos(par_theta) + orb['VZ'] * np.sin(par_theta))
	par_phi = np.transpose((orb['PHI'] + 2 * np.pi) % (2 * np.pi))
	par_vtheta = np.transpose(orb['VR'] * np.sin(par_theta) - orb['VZ'] * np.cos(par_theta))
	par_vphi = np.transpose(orb['VPHI'])
	
	if psimode is None:
		psimin = np.mean(rminor)
	else:
		psimin = psimode
	
	# here we solve for the exact resonant frequency based on orbit
	orbit_freq = abs(nnum * orb['CIQ']['W_PHI'] - mnum * orb['CIQ']['W_THETA']) / 2. / np.pi
	if fepm is None:
		fepm = [orbit_freq]
	
	power = t * 0.
	vdote_arr = np.zeros((len(psimin), len(fepm), len(phi0)))
	trend = np.copy(vdote_arr)
	
	last = 0
	npts = len(psimin) * len(fepm) * len(phi0)
	for rr in range(len(psimin)):
		if verbose:
			print(len(psimin) - rr)
		for ff in range(len(fepm)):
			for oo in range(len(phi0)):
				perc = (oo + ff * len(phi0) + len(fepm) * len(phi0) * rr) * 100. / (len(psimin) * len(fepm) * len(phi0))
				if perc - last > 1:
					print(f'{np.floor(perc)}% orbit done')
					last = np.floor(perc)
				
				psi_min_v = psimin[rr]
				fepm_v = fepm[ff]
				
				# todo: figure out psi version of this- don't have rhov, maybe interp2d?
				b_theta = np.interp(out['bp'], out['rhov'], psimin[rr])
				b_phi = np.interp(out['bphi'], out['rhov'], psimin[rr])
				q = np.interp(out['q'], out['rhov'], psimin[rr])
				
				phi = np.linspace(0, 2 * np.pi, endpoint=True)  # findgen(1001)/1000.*2*!pi
				theta = np.copy(phi)  # findgen(1001)/1000.*2*!pi
				
				amp0 = 1.
				sigma = .01  # guassian radial spread of mode
				
				# IF 0 THEN BEGIN ;as an aside, plot br contours
				# 	 rp = (findgen(1001)/1000.-.5)*1.04+1.5
				# 	 zp = (findgen(1001)/1000.-.5)*1.04
				#
				# 	 epm_pol = make_array(n_elements(rp),n_elements(zp))
				# 	 FOR r=0,n_elements(rp)-1 DO BEGIN
				# 	 FOR z=0,n_elements(zp)-1 DO $
				# 		epm_pol[r,z] = amp0*exp(-(sqrt((rp[r]-1.5)^2+zp[z]^2) - psimin) ^ 2 $
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
				
				field_line = (phi / q) % (2 * np.pi)
				
				Rphi = np.linspace(0, 2 * np.pi * r_mag)  # findgen(1001)/1000. * 2*!pi*r_maj
				rtheta = np.linspace(0, 2 * np.pi * psi_min_v)  # findgen(1001)/1000. * 2*!pi*r_min_v
				
				omega = 2 * np.pi * fepm_v
				delta_phi = omega * t + phi0[oo]
				line_pitch = b_theta / b_phi
				delta_theta = delta_phi * line_pitch
				
				# have delta_theta,delta_phi of mode as functions of time
				# as well as par_theta,par_phi location of particle as func of time
				
				amp = amp0 * np.exp(-(rminor - psimin[rr]) ^ 2 / 2. / sigma ^ 2)
				
				epm = amp * np.sin(nnum * par_phi + delta_phi) * np.sin(mnum * par_theta + delta_theta)
				
				# if 0:  # compute epm amp at particle location as func of time
				# 	power = epm*par_vr
				# 	vdote[rr,ff,oo] = int_tabulated(t,epm*par_vr)
				if 1:  # compute simplified vdote using curl (pg. 136)
					etheta = -1. / fepm[ff] * amp * np.cos(nnum * par_phi + delta_phi) * np.sin(
						mnum * par_theta + delta_theta)
					ephi = 1. / fepm[ff] * amp * np.sin(nnum * par_phi + delta_phi) * np.cos(
						mnum * par_theta + delta_theta)
					power = etheta * par_vtheta + ephi * par_vphi
					vdote_arr[rr, ff, oo] = integrate.simpson(power, t)
					if finescan:
						int_t = np.linspace(t[0], t[-1], endpoint=True, num=1001)
					else:
						int_t = np.linspace(t[0], t[-1], endpoint=True, num=1001)
					int_power = int_t * 0.
					for tt in range(len(int_t)):
						ii = max(
							np.where(t <= int_t[tt]))  # todo: this'll probably throw errors I always use np.where wrong
						int_power[tt] = integrate.simpson(power[ii], t[ii])
					fit = np.polyfit(int_t, int_power, 1)
					trend[rr, ff, oo] = fit[
						0]  # track slope of total power transfer (ORIGINAL USED abs(fit[0])- don't want that! want to see if power is going to particle or to mode)
		
		# end of offset array
	# end of fepm array
	# end of psimin array
	
	# todo: these steps to find index in psi, freq, phi probably won't work
	imax = np.where(trend == max(trend))
	omax = imax[0] / (len(psimin) * len(fepm))
	imax2 = imax[0] % (len(psimin) * len(fepm))
	fmax = imax2 / len(psimin)
	pmax = imax2 % len(psimin)
	
	if not singular:
		optimum = vdote(mnum, nnum, fepm[fmax], psimode[pmax], phi0[omax], verbose=False, plot=False, recon=recon,
		                orbit=orbit, nosave=True, replot=False, finescan=False)
		t = optimum.t
		power = optimum.power
		int_t = optimum.int_t
		int_power = optimum.int_power
	
	resonance = {'vdote': vdote_arr, 'psimode': psimin, 'fepm': fepm, 'phi0': phi0, 'fmax': fmax, 'pmax': pmax,
	             'omax': omax, 't': t, 'power': power, 'int_t': int_t, 'int_power': int_power, 'trend': trend}
	
	# PLOTTING
	vdote_arr = vdote_arr * 1.e5  # arbitrary scale anyway
	if plot and len(fepm) > 1 and len(psimin) > 1:
		fig, ax = plt.subplots()
		nl = 50
		lev = np.linspace(min(trend), max(trend), num=nl, endpoint=True)
		ax.contourf(trend[:, :, omax], psimin * 100., fepm / 1000.)
		ax.set_xlabel('r_mode (cm)')
		ax.set_ylabe('f (kHz')
		ax.set_zlabel('vdote (arb)')
		ax.plot(psimin[pmax] * 100., fepm[fmax] / 1000., 'o')
		if verbose:
			print(
				f'resonant at:\nr_mode = {psimin[pmax] * 100.}cm\nfepm = {fepm[fmax] / 1000.}kHz\nphi_0 = {phi0[omax]}')
	# !p.multi=[0,2,1]
	# yt = 6
	# ytn = replicate(' ',yt+1)
	# plotval = 1.1*[transpose(vdote[rmax,*,omax]),vdote[*,fmax,omax]]
	# yr = [min(plotval),max(plotval)]
	# plot,fepm/1000.,vdote[rmax,*,omax],psym=4,xtit='fepm (kHz)',yr=yr,/yst,$
	# 	 ytit='vdote (arb)',charsize=cs,yticks=yt,xmargin=[10,-3]
	# hline,0
	# plots,fepm[fmax]/1000.,vdote[rmax,fmax,omax],psym=4,color=!red,thick=2
	# billegend,['r_mode = ' + trim(psimin[rmax] * 100.)],/ right, box= 0, chars=cs
	# ;vdote vs r_mode
	# plot, psimin * 100., vdote[*, fmax, omax],psym= 4, xtit= 'r_mode (cm)', yr= yr,/ yst, $
	# 	 chars=cs,yticks=yt,ytickn=ytn,xmargin=[4,3]
	# hline,0
	# plots, psimin[rmax] * 100., vdote[rmax, fmax, omax], psym= 4, color=!red, thick=2
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
	# 	 print,'calculated rmode = '+trim(mean(rminor)*100.)+'cm'
	# 	 print,'calculated orbit freq = '+trim(orbit_freq/1000.)+'kHz'
	# 	ENDIF
	# ENDIF ;end of plotting
	
	return resonance


if __name__ == '__main__':
	pass
