
FUNCTION vdote, stop=stop, fitfile=fitfile,orbit=orbit, $
                mnum=mnum,nnum=nnum,nosave=nosave,fine=fine, $
                fepm=fepm,replot=replot,win=win,plot=plot,rmode=rmode, $
                verbose=verbose,phi0=phi0

  IF n_elements(fepm) EQ 1 AND n_elements(rmode) EQ 1 AND $
     n_elements(phi0) EQ 1 THEN singular = 1 ELSE singular = 0

  IF keyword_set(replot) THEN plot=1
  IF NOT keyword_set(win) THEN win = 0
  IF NOT keyword_set(replot) THEN BEGIN
       nf = 20.
       nr = 20.
       no = 1.
       IF NOT keyword_set(rmode) THEN rmode = (findgen(nr+1)/nr)*.25
       
       if not keyword_set(mnum) then m=1 ELSE m = mnum
       if not keyword_set(nnum) then n=5 ELSE n = nnum
       IF NOT keyword_set(fepm) THEN fepm = (findgen(nf+1)/nf-.5)*80.e3 + 100.e3
       IF NOT keyword_set(phi0) THEN phi0 = (findgen(no+1)/no)*2.*!pi

if not keyword_set(fitfile) THEN $
   fitfile = '/home/anderson/datasets/anderson/nbi_heating/bill_300kA/' + $
             'cinde_300kA_19.fit' ;Jay used these for qfi plot
restore, fitfile

if not keyword_set(orbit) then begin 
  ;get ion orbit: 
  dorb = findfile('/home/capecchi/resonance/orbits/d_phi*.sav')
  orbit = dorb(0)     
 ENDIF
restore,orbit 
t = transpose(data.t_orb) ;default orbit has 20,000 pts from 0-66.4microsec

;COMPUTE PARTICLE vr ARRAY (notebook2 pg 135)
r_maj = mean(data.r) ;average major radius of orbit (good est of shaf shift)
Rmaj_par = sqrt(data.x^2 + data.y^2)
rminor = sqrt( (Rmaj_par-r_maj)^2 + data.z^2)
par_theta = atan(data.z,Rmaj_par-r_maj)
par_theta = transpose((par_theta + 2*!pi) MOD (2*!pi))
par_vr = transpose(data.vr*cos(par_theta) + data.vz*sin(par_theta))
par_phi = transpose((data.phi + 10*2*!pi) MOD (2*!pi))
par_vtheta = transpose(data.vr*sin(par_theta) - data.vz*cos(par_theta))
par_vphi = transpose(data.vphi)

IF NOT keyword_set(rmode) THEN r_min = mean(rminor) ELSE r_min = rmode

;here we solve for the exact resonant frequency based on orbit
orbit_freq = abs( n* data.ciq.w_phi - m*data.ciq.w_theta)/2./!pi 
if not keyword_set(fepm) then fepm = orbit_freq
;if fepm lt 0 then fepm=-fepm 
;ENDIF

power = t*0.
vdote = make_array(n_elements(r_min),n_elements(fepm),n_elements(phi0))
trend = vdote

last = 0
npts = n_elements(r_min)*n_elements(fepm)*n_elements(phi0)
FOR rr=0,n_elements(r_min)-1 DO BEGIN
     IF keyword_set(verbose) THEN print,n_elements(r_min)-rr
FOR ff=0,n_elements(fepm)-1 DO BEGIN
FOR oo=0,n_elements(phi0)-1 DO BEGIN
     perc = (oo+ff*n_elements(phi0)+n_elements(fepm)*n_elements(phi0)*rr)*100./(n_elements(r_min)*n_elements(fepm)*n_elements(phi0))
     IF perc-last GE 1 THEN BEGIN
          print,trim(fix(perc))+'% orbit done'
          last = fix(perc)
         ENDIF

r_min_v = r_min[rr]
fepm_v = fepm[ff]

;iloc = closest(out.rhov,r_min_v)
;b_theta = out.bp(iloc(0)) 
;b_phi = out.bphi(iloc(0))
;q = out.q(iloc(0))
b_theta = interpol(out.bp,out.rhov,r_min[rr])
b_phi = interpol(out.bphi,out.rhov,r_min[rr])
q = interpol(out.q,out.rhov,r_min[rr])

phi = findgen(1001)/1000.*2*!pi
theta = findgen(1001)/1000.*2*!pi

amp0 = 1.
sigma = .01 ;guassian radial spread of epm

IF 0 THEN BEGIN ;as an aside, plot br contours
     rp = (findgen(1001)/1000.-.5)*1.04+1.5
     zp = (findgen(1001)/1000.-.5)*1.04
     
     epm_pol = make_array(n_elements(rp),n_elements(zp))
     FOR r=0,n_elements(rp)-1 DO BEGIN
     FOR z=0,n_elements(zp)-1 DO $
        epm_pol[r,z] = amp0*exp(-(sqrt((rp[r]-1.5)^2+zp[z]^2)-r_min)^2 $
                       /2./sigma^2)*sin(m*atan(zp[z],rp[r]-1.5)+!pi/3.)
    ENDFOR
     th = findgen(1001)/1000.*!pi
     epm_tor = make_array(n_elements(th),n_elements(phi))
     FOR p=0,n_elements(phi)-1 DO BEGIN
          FOR t=0,n_elements(th)-1 DO BEGIN
               epm_tor[t,p] = amp0*sin(m*th[t])*sin(n*phi[p])
              ENDFOR
         ENDFOR
     
     putwindow,5
     !p.multi=[0,2,1]
     nl=50
     lev=findgen(nl)/(nl-1)*(max(epm_pol)-min(epm_pol))+min(epm_pol)
     contour,epm_pol,rp,zp,xtit='radius (m)',ytit='z (m)', $
             lev=lev,/iso,/fill
     contour,epm_tor,th/!pi,phi/!pi,levels=findgen(25)/24.*(max(epm_tor)- $
             min(epm_tor))+min(epm_tor),/fill,xr=[0,1.],yr=[0,2.], $
             xtit='theta/pi',ytit='phi/pi'

ENDIF ;end of plot br contours

field_line = (phi / q ) mod (2*!pi)

Rphi = findgen(1001)/1000. * 2*!pi*r_maj
rtheta = findgen(1001)/1000. * 2*!pi*r_min_v

;offset = !pi
omega = 2*!pi*fepm_v
delta_phi = omega*t+ phi0[oo]
line_pitch = b_theta/b_phi
delta_theta = delta_phi*line_pitch

;have delta_theta,delta_phi of mode as functions of time
;as well as par_theta,par_phi location of particle as func of time

amp = amp0*exp(-(rminor-r_min[rr])^2/2./sigma^2)

epm = amp*sin(n*par_phi+delta_phi)*sin(m*par_theta + delta_theta)

IF 0 THEN BEGIN ;compute epm amp at particle location as func of time
power = epm*par_vr
vdote[rr,ff,oo] = int_tabulated(t,epm*par_vr)
ENDIF
IF 1 THEN BEGIN ;compute simplified vdote using curl (pg. 136)
etheta = -1./fepm[ff]*amp*cos(n*par_phi+delta_phi)*sin(m*par_theta+delta_theta)
ephi = 1./fepm[ff]*amp*sin(n*par_phi+delta_phi)*cos(m*par_theta+delta_theta)
power = etheta*par_vtheta + ephi*par_vphi
vdote[rr,ff,oo] = int_tabulated(t,etheta*par_vtheta + ephi*par_vphi)
IF NOT keyword_set(fine) THEN $
   int_t = findgen(101)/100.*(t[-1]-t[0])+t[0] ELSE $
      int_t = findgen(1001)/1000.*(t[-1]-t[0])+t[0]
int_power = int_t*0.

FOR tt=1,n_elements(int_t)-1 DO BEGIN
     ii = where(t LE int_t[tt])
     int_power[tt] = int_tabulated(t[ii],power[ii])
    ENDFOR
fit = linfit(int_t,int_power)
trend[rr,ff,oo] = abs(fit[1]) ;track slope of total power transfer
ENDIF 

ENDFOR                          ;end of offset array
ENDFOR                          ;end of fepm array
ENDFOR                          ;end of r_min array

; imax = where(vdote EQ max(vdote))
 imax = where(trend EQ max(trend))
 omax = imax[0]/(n_elements(r_min)*n_elements(fepm))
 imax2 = imax[0] MOD (n_elements(r_min)*n_elements(fepm))
 fmax = imax2/n_elements(r_min)
 rmax = imax2 MOD n_elements(r_min)

 IF NOT singular THEN BEGIN
      optimum = vdote(fitfile=fitfile,orbit=orbit,mnum=mnum,nnum=nnum, $
                      fepm=fepm[fmax],rmode=rmode[rmax],phi0=phi0[omax], $
                      nosave=nosave)
      t = optimum.t
      power = optimum.power
      int_t = optimum.int_t
      int_power = optimum.int_power
     ENDIF 

 IF NOT keyword_set(nosave) THEN $
    save,filename='/home/capecchi/resonance/resonance_data.dat'
ENDIF ELSE restore,'/home/capecchi/resonance/resonance_data.dat' ;end of replot


;PLOTTING
vdote = vdote*1.e5         ;arbitrary scale anyway
cs=1.5
IF keyword_set(plot) THEN BEGIN
     IF n_elements(fepm) GT 1 AND n_elements(r_min) GT 1 THEN BEGIN
 putwindow,win
 nl=50
 IF 0 THEN BEGIN ;vdote
 lev=findgen(nl)/(nl-1)*(max(vdote)-min(vdote))+min(vdote)
 contour_bar,vdote[*,*,omax],r_min*100.,fepm/1000.,xtit='r_mode (cm)', $
         ytit='fepm (kHz)',ztit='vdote (arb)',chars=cs,lev=lev,/fill
ENDIF ELSE BEGIN
 lev=findgen(nl)/(nl-1)*(max(trend)-min(trend))+min(trend)
 contour_bar,trend[*,*,omax],r_min*100.,fepm/1000.,xtit='r_mode (cm)', $
         ytit='fepm (kHz)',ztit='vdote (arb)',chars=cs,lev=lev,/fill
ENDELSE
 plots,r_min[rmax]*100.,fepm[fmax]/1000.,psym=4,color=0,thick=2
 IF keyword_set(verbose) THEN BEGIN 
      print,'resonant at:'
      print,'r_mode = '+trim(r_min[rmax]*100.)+'cm'
      print,'fepm = '+trim(fepm[fmax]/1000.)+'kHz'
      print,'phi_0 = '+trim(phi0[omax])
     ENDIF ;end of verbose
ENDIF

putwindow,1
!p.multi=[0,2,1]
yt = 6
ytn = replicate(' ',yt+1)
plotval = 1.1*[transpose(vdote[rmax,*,omax]),vdote[*,fmax,omax]]
yr = [min(plotval),max(plotval)]
plot,fepm/1000.,vdote[rmax,*,omax],psym=4,xtit='fepm (kHz)',yr=yr,/yst,$
     ytit='vdote (arb)',charsize=cs,yticks=yt,xmargin=[10,-3]
hline,0
plots,fepm[fmax]/1000.,vdote[rmax,fmax,omax],psym=4,color=!red,thick=2
billegend,['r_mode = '+trim(r_min[rmax]*100.)],/right,box=0,chars=cs
;vdote vs r_mode
plot,r_min*100.,vdote[*,fmax,omax],psym=4,xtit='r_mode (cm)',yr=yr,/yst, $
     chars=cs,yticks=yt,ytickn=ytn,xmargin=[4,3]
hline,0
plots,r_min[rmax]*100.,vdote[rmax,fmax,omax],psym=4,color=!red,thick=2
billegend,['fepm = '+trim(fepm[fmax]/1000.)],/right,box=0,chars=cs

putwindow,2
!p.multi=[0,1,2]
xt = 8
xtn = replicate(' ',xt+1)
t = t*1.e6 ;convert to microseconds
int_t = int_t*1.e6
plot,t,power,ytit='power(t) (arb)', $
     chars=cs,xticks=xt,xtickn=xtn,ymar=[2,2]
int_power = int_power*1.e5 ;arbitrary anyway
plot,int_t,int_power,ymar=[4,-1],xticks=xt,chars=cs, $
     ytit=textoidl('int(P)'),xtit=textoidl('time (\mu s)')

putwindow,4
!p.multi=[0,2,1]
plot,data.r[0:10000],data.z[0:10000],/iso,xtit='R (m)',ytit='Z (m)'
plot,data.x[0:10000],data.y[0:10000],/iso,xtit='X (m)',ytit='Y (m)'

IF keyword_set(verbose) THEN BEGIN 
     print,'calculated rmode = '+trim(mean(rminor)*100.)+'cm'
     print,'calculated orbit freq = '+trim(orbit_freq/1000.)+'kHz'
    ENDIF 
ENDIF ;end of plotting

;IF singular THEN $

resonance = {vdote:vdote,rmode:r_min,fepm:fepm,phi0:phi0,fmax:fmax,rmax:rmax, $
             omax:omax,t:t,power:power,int_t:int_t,int_power:int_power, $
             trend:trend}
IF NOT keyword_set(nosave) THEN $
   save,filename=orbit,input,data,resonance

IF keyword_set(stop) THEN stop
return,resonance
END
