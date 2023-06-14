;--------------------------------------------------------------
; Interpolates between known values for the charge exchange area
function xsection,ea;vplasma
common constants_cgs, mp,kb,kev_erg,AU_cm
;modify later to use Lindsay and Stebbings (2005) cross sections formula  
;Input:
; Ea = total energy in keV of H+ (soon to be ENA)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Old Prested cross-section technique
;
;Prob=[[0.0, 30.50],[0.3, 20.0],[7.0, 10.0],[10.01, 8.0],$
;[20.0, 4.0],[30.0, 3.0],[40.0, 2.0],[50.0, 1.0],$
;[60.0, 0.5],[70.0, 0.4],[80.0, 0.3],[90.0, 0.2],$
;[100.0, 0.1],[200.0, 0.01],[300.0, 0.001],[400.0, 3.0*10.0^(-4.0)],$
;[500.0, 0.0001],[600.0, 7.0*10.0^(-5.0)],[700.0, 2.0*10.0^(-5.0)]];   Prob[0]=energy keV, Prob[1]= sigma cm^-2

; sigma may be off by 13.6 eV due to losses in charge exchange, negliable
;if ea lt prob[0,0] then ea = prob[0,0]
;if ea gt prob[0,18] then ea = prob[0,18]

;for i=0, 17 do begin
;    if ea ge prob[0,i] and ea le prob[0,i+1] then begin
;        sigma=prob[1,i]+((ea-prob[0,i])/(prob[0,i+1]-prob[0,i]))*(prob[1,i+1]-prob[1,i])
;    endif
;endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if ea ne ea or ea lt 0 then stop

; Lindsay & Stebbings (2005) cross section formula - added by MZK
a1=4.15
a2=0.531
a3=67.3

if ea ne ea or ea lt 0 then stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Lindsay & Stebbings (2005) cross-section from relative velocity
; inspired by Bzowski (2008), Appendix A.1
;en=0.5*mp*vplasma^2.
;en=en/kev_erg
;sig1=(a1-a2*alog(en))^2.
;sig2=(1-exp(-a3/en))^(4.5)
;sigma=sig1*sig2
;sigma=sigma*10D^(-16D)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Lindsay & Stebbings (2005) cross-section from ENA velocity
sig1=(a1-a2*alog(ea))^2.
sig2=(1-exp(-a3/ea))^(4.5)
sigma=sig1*sig2
sigma=sigma*10D^(-16D)

if ea ne ea or ea lt 0 then stop

;sigma=1.29*10D^(-17D) ; for Helium @ 100 keV
;sigma=2.33*10D^(-16D) ; for Helium @ 10 keV

if sigma ne sigma then stop
return, sigma
end
;----------------------------------------------------------------
;Finds the unitless probability of an ena of some energy to reach the earth without getting re-ionized
;Can make more complicated
function SurvivalProb,Ea,r,theta_value,vplasma,sigma,np,bet
;Inputs:
; Ea = Energy of ENA in keV 
; R = Radial distance in AU
common constants_cgs, mp,kb,kev_erg,AU_cm
common grid_params, dr, dtheta ; dtheta by MZK

R1=AU_cm   ;earth's radial distance in cm
Va=energy_velocity(Ea)    ;cm/s
Beta0 = 6.0*10.0^(-7.0)  ;~ ionization rate in s^-1

;S=exp(-1.0*Beta0*r1*(1.-(1./r))*va^(-1.0))
;S=exp(-1.0*Beta0*r1*(1.-(1./r))*!pi/2.*va^(-1.0)); - IBEX path from Heerikhuisen et al. (2008) - MZK
S=exp(-1.0*bet*r1/va) ; integrated survival probability - MZK

;-----------------------------------------------
;For Bzowski survivability profile
;theta_value = theta_value - 90. 
;For 0.1 keV
;From figure 17 upper panel, in M. Bzowski (2008)
;theta_arr=[-90., -75., -65., -55., -50., -45., -40., -35., -30., -25., -5., 5., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90.]
;Sarr    = [ .43, .425, .415,  .4,  .385,  .37,  .35,  .33, .325,  .32, .32,.32,.325, .33, .35, .38, .42, .45, .47, .48,.495,  .5,.505, .51]

;for 1 kEv
;Extrapolated from Fig 15 (Right Lower) and using Gruntman survivalbility for 0 degrees (.82)
;theta_arr    =    [-90., -75., -65., -55., -50., -40., -30., 15., 25., 30., 35., 40., 90.]
;Sarr         =0.82*[1.12, 1.11,1.105,  1.1, 1.08, 1.05,   1.,  1.,1.05, 1.1,1.13,1.15, 1.2]


;for i=0, 11 do begin;22 do begin ;22 tot for .1 , 11 for 1. keV
;    if theta_value ge theta_arr[i] and theta_value le theta_arr[i+1] then begin
;        S=Sarr[i]+((theta_value-theta_arr[i])/(theta_arr[i+1]-theta_arr[i]))*(Sarr[i+1]-Sarr[i])
;    endif
;endfor

;------------------------------------------------
if S ne S then stop
return, S
end

;----------------------------------------------------------------
;Multiplies together all the constant factors that go into finding
;flux (current)

function fluxconstants,Ea,r,nn,theta_value,vplasma,np,bet
common include_survival, flag_survival
;Inputs: 
; nn = neutral number density in cm^-3
; r = distance from sun in AU
; Ea = Energy of ENA in keV
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen
common constants_cgs, mp,kb,kev_erg,AU_cm

if ea ne ea or r ne r or nn ne nn or theta_value ne theta_value or vplasma ne vplasma or bet ne bet then stop
if ea lt 0 then stop

L_los = dr
;L_los=1.0
L_loscm=L_los*AU_cm; convert from au to cm
sigma=xsection(ea)
S = SurvivalProb(Ea,r,theta_value,vplasma,sigma,np,bet)
;S = SurvivalProb(Ea,r,theta_value,vplasma,sigma,np)

if flag_survival eq 1 then fluxconstant=sigma*L_Loscm*nn*S else fluxconstant=sigma*L_Loscm*nn  ; 

if fluxconstant ne fluxconstant then stop
return,fluxconstant
end
                           
