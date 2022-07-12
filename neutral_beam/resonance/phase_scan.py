import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from neutral_beam.resonance.resonance_tools import convert_sav_to_dict
import pickle
from neutral_beam.resonance.vdote import vdote

if os.path.exists('Z:/users/wcapecch/'):
	direc = 'Z:/users/wcapecch/'
	proj_direc = 'Z:/PycharmProjects/LTXb-py/'
elif os.path.exists('//samba/wcapecch/'):
	direc = '//samba/wcapecch/'
	proj_direc = '//samba/wcapecch/PycharmProjects/LTXb-py/'


def phase_scan_ltx(mnum=1, nnum=2, fepm=None, psimode=None, phi0=None, verbose=False):  # , parallel=parallel,
	# ps=ps, click=click, newdata=newdata, contours=contours, thesis=thesis, reclick=reclick)
	
	nf = 100
	nr = 100
	no = 50
	
	if psimode is None:
		psimode = np.linspace(0., .35, num=nr)
	if fepm is None:
		fepm = np.linspace(25.e3, 175.e3, num=nf)
	if phi0 is None:
		phi0 = 0.
	
	orbit_dir = f'{proj_direc}neutral_beam/resonance/orbits/'
	orbits = glob.glob(f'{orbit_dir}*orbit.sav')
	for orbit in orbits:
		if not os.path.isfile(f'{orbit[:-4]}.p'):  # pickled version doesn't exist
			convert_sav_to_dict(orbit)
	orbits = glob.glob(f'{orbit_dir}*orbit.p')  # list of pickled orbits
	
	r_res = []  # array of computed resonance radii
	f_res = []  # same for freqs of epm
	o_res = []  # same for phi0 of epm
	r_orb = []  # array of orbit computed radii
	f_orb = np.zeros((3, len(orbits)))  # orbit predicted fepm for l=-1,0,1
	ioniz0 = np.zeros((4, len(orbits)))  # x,y,R,z of ionization pt for each orbit
	gc0 = np.zeros((4, len(orbits)))  # x,y,R,z of guiding center pt for each orbit at ionization
	
	for o in range(len(orbits)):
		orbit_dat = pickle.load(open(orbits[o], 'rb'))
		orb = orbit_dat['orb']
		ioniz0[:, o] = [orb['X'][0], orb['Y'][0], orb['R'][0], orb['Z'][0]]
		# TODO compute GC x,y,R,z and update gc0 array- DO I NEED GC DATA? maybe for plotting initial GC at ionization?
		print(orbits[o])
		if 'resonance' not in orbit_dat.keys():
			print(f'{len(orbits) - o} orbits remaining')
			recon = orbit_dat['orb']['RECONST_FILE'].decode()
			resonance = vdote(mnum, nnum, fepm, psimode, phi0, verbose=verbose, plot=False, recon=recon,
			                  orbit=orbits[o])
			orbit_dat['resonance'] = resonance
			pickle.dump(orbit_dat, open(orbits[o], 'wb'))  # add resonance data to orbit_dat file
		r_mag = orb['OUT']['R0']  # [m]
		Rmaj_par = np.sqrt(orb['X'] ** 2 + orb['Y'] ** 2)
		rminor = np.sqrt((Rmaj_par - r_mag) ** 2 + orb['Z'] ** 2)
		of_l0 = (nnum * orb['CIQ']['W_PHI'] - mnum * orb['CIQ']['W_THETA']) / 2. / np.pi
		of_lm1 = (nnum * orb['CIQ']['W_PHI'] - (mnum - 1) * orb['CIQ']['W_THETA']) / 2. / np.pi
		of_lp1 = (nnum * orb['CIQ']['W_PHI'] - (mnum + 1) * orb['CIQ']['W_THETA']) / 2. / np.pi
		orbit_r = np.mean(rminor)
		
		o_res.append(resonance['phi0'][resonance['omax']])
		r_orb.append(orbit_r)
		f_orb[:, o] = [of_lm1, of_l0, of_lp1]
		r_res.append(resonance['rmode'][resonance['rmax']])
		f_res.append(resonance['fepm'][resonance['fmax']])
	# 	IF keyword_set(contours) OR keyword_set(reclick) THEN BEGIN
	# 		IF NOT keyword_set(click) AND NOT keyword_set(ps) THEN putwindow,o
	# 		  IF keyword_set(reclick) THEN putwindow,0
	# 		  nl=50
	# 		  IF 0 THEN BEGIN       ;vdote
	# ; vdote = resonance.vdote[*,*,resonance.omax]
	# ; lev=findgen(nl)/(nl-1)*(max(vdote));-min(vdote))+min(vdote)
	# ; contour_bar,resonance.vdote[*,*,resonance.omax],resonance.rmode*100., $
	# ;             resonance.fepm/1000.,xtit='r_mode (cm)',ytit='fepm (kHz)', $
	# ;             tit='vdote (arb)',chars=cs,lev=lev,/fill
	# ; plots,resonance.rmode[resonance.rmax]*100., $
	# ;       resonance.fepm[resonance.fmax]/1000.,psym=4,color=0,thick=2
	# 			  ENDIF ELSE BEGIN
	# 				   !p.multi=[0,1,2]
	# 				   trend = resonance.trend[*,*,resonance.omax]
	# 				   lev=findgen(nl)/(nl-1)*(max(trend)) ;-min(trend))+min(trend)
	# 				   contour,trend,resonance.rmode,ymar=[-2,2], $
	# 						   resonance.fepm,xtit='r_mode (m)',chars=cs,lev=lev, $
	# 						   /fill,ytit='fepm (Hz)',yr=[25,175]*1000.,/ysty,xr=[0,35]/100., $
	# 						   /xsty,tit='Power Transfer Trend (arb)'
	# 				   IF NOT keyword_set(reclick) THEN BEGIN
	# 						imax = where(trend EQ max(trend))
	# 						fmax = imax[0]/n_elements(resonance.rmode)
	# 						rmax = imax[0] MOD n_elements(resonance.rmode)
	# 						int_t = resonance.int_t
	# 						int_power = resonance.int_power
	# 					   ENDIF ELSE BEGIN
	# 							cursor,rval,fval,/data,/down
	# 							rmax = max(where(resonance.rmode*100. LE rval))
	# 							fmax = max(where(resonance.fepm/1000. LE fval))
	# 							clickres = vdote(mnum=mnum,nnum=nnum,fepm=fval, $
	# 											 rmode=rval,recon=recon,/nosave, $
	# 											 orbit=orbits[o],phi0=phi0)
	# 							int_t = clickres.int_t
	# 							int_power = clickres.int_power
	# 							save,filename=orbits[o],input,data,resonance,clickres
	# 						   ENDELSE
	# 					   plots,resonance.rmode[rmax], $
	# 							 resonance.fepm[fmax],psym=4,color=0,thick=2
	# 					   plot,int_t,int_power,ymar=[4,6]
	#
	# ;IF NOT keyword_set(ps) THEN putwindow,o+2
	# ; contour_bar,trend,resonance.rmode*100.,resonance.fepm/1000.
	#
	# ;abs(resonance.trend[*,*,resonance.omax-50]),resonance.rmode*100.,resonance.fepm/1000.
	#
	# 					  ENDELSE   ;end
	# 			 ENDIF              ;end of contours
	#
	# 	 IF NOT keyword_set(click) THEN BEGIN
	# 		  r_res = [r_res,resonance.rmode[resonance.rmax]]
	# 		  f_res = [f_res,resonance.fepm[resonance.fmax]]
	# 		 ENDIF ELSE BEGIN
	# 			  r_res = [r_res,clickres.rmode]
	# 			  f_res = [f_res,clickres.fepm]
	# 			 ENDELSE
	# 		ENDFOR ; end of orbits
	
	# PLOTTING
	# ;restore,'/home/capecchi/resonance/resonance_plot.dat';bc,indexes,phi0
	#    IF keyword_set(thesis) THEN BEGIN
	# 		fthesis = '/home/capecchi/thesis_plots/phase_scan'
	# 		thesis_plot_start,file=fthesis,/not_mstfit
	# 	   ENDIF ELSE putwindow,0 ;fepm * f_orbit vs rmode
	# cs=1.25
	# clrs = fix(findgen(n_elements(r_res))/(n_elements(r_res)-1)*254)
	
	# topx = (findgen(2)-.5)*1.*2.02*2
	# topy = (findgen(2)-.5)*1.2*2.02*2
	# xr = (findgen(2)-.5)*1.*.52*2+1.5
	# yr = (findgen(2)-.5)*1.2*.52*2
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
	# !p.multi=[0,4,2]
	theta = np.linspace(0, 2 * np.pi, endpoint=True)
	# rr = 1.5+.52*cos(theta)
	# zz = .52*sin(theta)
	# xout = 2.02*cos(theta)
	# xin = .98*cos(theta)
	# yout = 2.02*sin(theta)
	# yin = .98*sin(theta)
	# choose an orbit to display as example
	chosen = 3
	pchosen = 2
	
	# draw_ltx top
	# ax1.plot(min(orb['OUT']['RVES'])*np.cos(theta),min(orb['OUT']['RVES'])*np.sin(theta),'k-')
	# ax1.plot(max(orb['OUT']['RVES'])*np.cos(theta),max(orb['OUT']['RVES'])*np.sin(theta),'k-')
	# oplot,bc.x,bc.y,linestyle=2
	# FOR p=0,n_elements(r_res)-1 DO $
	#    plots,bc.x[iorder[p]],bc.y[iorder[p]],psym=4,color=clrs[p],thick=1.5
	# plots,bc.x[iorder[chosen]],bc.y[iorder[chosen]],psym=pchosen, $
	# 	  color=clrs[chosen],thick=1.5,symsize=1.5
	
	# draw ltx cross-section
	ax1.plot(orb['OUT']['RVES'], orb['OUT']['ZVES'], 'k-')
	
	# plot,rr,zz,/iso,xr=xr,yr=yr,xstyle=5,ystyle=5,xmar=[4,3]
	# oplot,bc.r,bc.z,linestyle=2
	# FOR p=0,n_elements(r_res)-1 DO $
	#    plots,bc.r[iorder[p]],bc.z[iorder[p]],psym=4,color=clrs[p],thick=1.5
	# plots,bc.r[iorder[chosen]],bc.z[iorder[chosen]],psym=pchosen, $
	# 	  color=clrs[chosen],thick=1.5,symsize=1.5
	
	# !p.multi=[3,2,2]
	# restore,orbits[chosen] ;dunno which one this is
	
	# spectrogram for chosen orbit
	orbit_dat = pickle.load(open(orbits[chosen], 'rb'))
	trend = orbit_dat['resonance']['trend'][:, :, orbit_dat['resonance']['omax']]
	rmode = orbit_dat['resonance']['rmode'] * 100.
	fepm = orbit_dat['resonance']['fepm'] / 1000.
	nl = 50
	# lev=findgen(nl)/(nl-1)*(max(trend));-min(trend))+min(trend)
	ax2.contourf(trend, rmode, fepm)
	ax2.set_xlabel('r_mode (cm)')
	ax2.set_ylabel('f (kHz)')
	ax2.set_title('Power Transfer Trend (arb)')


# fepm vs rmode
# xr = [0.,25.]
# yr = [0,150.]
# clrs2 = clrs
# clrs2 = clrs*0
# iok = where(f_res LE 150000.)
# f_orb0 = transpose(f_orb[1,*])
#    plot,r_res[iok]*100.,f_res[iok]/1000.,/nodata, $
# 		xtit='rmode (cm)',ytit='fepm (kHz)', $
# 		chars=cs,xr=xr,yr=yr
#    FOR p=0,n_elements(iok)-1 DO $
# 	  plots,r_res[iok[p]]*100.,f_res[iok[p]]/1000.,psym=4,color=clrs[iok[p]]
#    plots,r_res[chosen]*100.,f_res[chosen]/1000., $
# 		 psym=pchosen,color=clrs[chosen],symsize=1.5
#
# ;power trend for chosen orbit
#    int_t = clickres.int_t
#    int_power = clickres.int_power*1.e6 ;arbitrary scale anyway
#    plot,int_t,int_power,xtit='time (s)',ytit='Power Transferred (arb)', $
# 		chars=cs
#
#    IF keyword_set(thesis) THEN BEGIN
# 		save,filename='/home/capecchi/state.dat'
# 		used = ['r_res','f_res','f_orb0','trend','fepm','rmode', $
# 				'int_t','int_power']
# 		thesis_plot_finish,file=fthesis,/noprompt,/csv,savethese=used
# 	   ENDIF
#
# ;print resonance parameters for chosen orbit;
# print,' '
# print,'rmode:: '+trim(r_res[chosen])
# print,'fepm:: '+trim(f_res[chosen])


if __name__ == '__main__':
	phase_scan_ltx()
	pass
