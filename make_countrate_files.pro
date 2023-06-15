;----------------------------------------------


;Finds the Energetic Neutral atom flux for specific place in the
;heliosphere and distribution
;has to be first to complie correctly

function calc_diff_flux,dist,Ea,r,Tp,Ur,Un,nn,np,theta,kappa0,bet ; kappa0 added by MZK
;-Inputs
;Energy of desired ENA (keV)
;Radius (AU)
; Line of sight (los integral is brokeninto reimann sum) (AU)
;temperature of the plasma (K)
;Bulk radial velocity of plasma (cm/s)
;bulk tangential speed of plasma (cm/s),
;neutral density (cm^-3)
;plasma density (cm^-3)
;-Output
;Current of ENAs'Kappa'+string(kappa,format='(f4.2)')
common kappa_value,kappa,kappaPUI,kappaRef,kappaISM,kappaOH
common constants_cgs, mp,kb,kev_erg, AU_cm
common distributions,dist_names
;put energy into solar wind frame and convert to velocity

if kappa0 le 1.50000 and dist eq 'shell' then dist='Maxwellian'
if kappa0 le 1.50000 and dist eq 'Kappa'+string(kappa,format='(f4.2)') then dist='Maxwellian'
if kappa0 gt 1.50000 and dist eq 'Kappa'+string(kappa,format='(f4.2)') then dist='Kappa'+string(kappa,format='(f4.2)')
v_plasma = plasma_frame(Ea,Ur,Un)


if dist eq 'Kappa'+string(kappa,format='(f4.2)') then f=double(kappa_fn(Tp,np,v_plasma,kappa0)) ; kappa0 added by MZK
if dist eq 'Maxwellian' then f=double(maxwellian_fn(Tp,np,v_plasma))
if dist eq 'shell' then f=double(shell_fn(np,v_plasma,Ur)); from Zank et al. (2010)

if ea ne ea or r ne r or nn ne nn or theta ne theta or v_plasma ne v_plasma or bet ne bet then stop

fluxconstant=double(fluxconstants(Ea,r,nn,theta,v_plasma,np,bet)) ;Sigma*Nn*S*survival probability

j=double((2.0*Ea*kev_erg*f*fluxconstant)/(mp^2.0)) ;s cm^-4 sr^-1 g^-1= 1/(cm^2*g/s^2) s^-1 sr^-1 cm^-2 (cgs)

j=j*kev_erg ;converts from erg^-1 to kev^-1

if j ne j then stop
dist=dist_names(0)

return, j
end
;---------------------------------------

;Finds the Energetic Neutral atom flux for specific place in the
;heliosphere and distribution
;has to be first to complie correctly

function calc_diff_flux_zank,dist,Ea,r,Ur,Un,nn,np,theta,del,bet,el,em,vshock ; kappa0 added by MZK
;-Inputs
;Energy of desired ENA (keV)
;Radius (AU)
; Line of sight (los integral is brokeninto reimann sum) (AU)
;temperature of the plasma (K)
;Bulk radial velocity of plasma (cm/s)
;bulk tangential speed of plasma (cm/s),
;neutral density (cm^-3)
;plasma density (cm^-3)
;-Output
;Current of ENAs'Kappa'+string(kappa,format='(f4.2)')
common kappa_value,kappa,kappaPUI,kappaRef,kappaISM,kappaOH
common constants_cgs, mp,kb,kev_erg, AU_cm
common distributions,dist_names
;put energy into solar wind frame and convert to velocity

dist='Zank'
v_plasma = plasma_frame(Ea,Ur,Un)-vshock ; accounting for shock frame in which distribution is calculated

f=double(zank_fn(v_plasma,np,del,el,em,Ea))

if ea ne ea or r ne r or nn ne nn or theta ne theta or v_plasma ne v_plasma or bet ne bet then stop

fluxconstant=double(fluxconstants(Ea,r,nn,theta,v_plasma,np,bet)) ;Sigma*Nn*S*survival probability

j=double((2.0*Ea*kev_erg*f*fluxconstant)/(mp^2.0)) ;s cm^-4 sr^-1 g^-1= 1/(cm^2*g/s^2) s^-1 sr^-1 cm^-2 (cgs)

j=j*kev_erg ;converts from erg^-1 to kev^-1

if j ne j then stop
dist=dist_names(0)

return, j
end

;---------------------------------------

pro allCountRate
;Creates count rate data files for all energy bands for a specific
;list of distribution functions and a specific list of modelers
;(Merav Opher or Jacob Heerkhuissen).
;Outputs:
; saves countrate maps to files

;can't imagine needing resolution smaller than 10 eV, could be adjusted
eres = 10.

common distributions,dist_names,num_dist
common model_names, models, num_models, pui_model,vasyliunas,ribbon,gdf,hinterp,ntr_frac,nref_frac,ttr_frac,tref_frac,colde,esw_frac,etr_frac,eref_frac,moscow,moscow_temp,giacalone
common sensor_params,hi,lo, wt9, wt10, wt11, wt12, wt13, ultra
common all_sensors,Eminarr,Emaxarr,Garr,Ecenarr
common energybands,num_ebands,arr_eband_indexes
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen,r_lower,r_upper
common grid_opher,gridopher
common include_survival, flag_survival, extinction

nr=gridopher.r.num

;makes one array out of hi and lo

Eminarr=[lo.E.min,hi.e.min,ultra.e.min]
Emaxarr=[lo.e.max,hi.e.max,ultra.e.max]
Garr=[lo.g,hi.g,ultra.g]
Ecenarr=[lo.e.central,hi.e.central,ultra.e.central]

lc=make_array(ntheta,nphi)
cool=make_array(ntheta,nphi,nr,/double)
vx=make_array(ntheta,nphi,nr,/double)
vy=make_array(ntheta,nphi,nr,/double)
vz=make_array(ntheta,nphi,nr,/double)
a_z=make_array(ntheta,nphi,nr,/double)
np_stream=make_array(ntheta,nphi,nr,/double)
tp_stream=make_array(ntheta,nphi,nr,/double)
vr_stream=make_array(ntheta,nphi,nr,/double)
vt_stream=make_array(ntheta,nphi,nr,/double)
ntr_frac0=make_array(ntheta,nphi,nr,/double)
nref_frac0=make_array(ntheta,nphi,nr,/double)
nen_frac0=make_array(ntheta,nphi,nr,/double)
etr_frac0=make_array(ntheta,nphi,nr,/double)
eref_frac0=make_array(ntheta,nphi,nr,/double)
een_frac0=make_array(ntheta,nphi,nr,/double)
em_frac0=make_array(ntheta,nphi,nr,/double)
kappa_ref0=make_array(ntheta,nphi,nr,/double)
vshock0=make_array(ntheta,nphi,nr,/double)
esw_frac0=make_array(ntheta,nphi,nr,/double)

if extinction eq 1 and giacalone eq 0 or vasyliunas eq 1 then read_stream, lc, cool, vx, vy, vz, a_z, np_stream, tp_stream, vr_stream, vt_stream ; GIACALONE
if giacalone eq 1 then read_stream_g, lc, cool, vx, vy, vz, ntr_frac0, nref_frac0, nen_frac0, etr_frac0, eref_frac0, een_frac0,em_frac0, kappa_ref0,vshock0,esw_frac0 ; GIACALONE 

;do for all modelers and all distributions
for i=0,num_models-1 do begin
    for j=0,num_dist-1 do begin
            for k=0.,num_ebands-1 do begin
                   ebandindex =arr_eband_indexes(k)
;                   countrate,ebandindex,eres,models(i),dist_names(j) ; original w/o cooling depth
                   countrate,ebandindex,eres,models(i),dist_names(j), lc, cool, vx, vy, vz, a_z, np_stream, tp_stream, vr_stream, vt_stream, $ 
                   ntr_frac0, nref_frac0, nen_frac0, etr_frac0, eref_frac0, een_frac0, em_frac0, kappa_ref0,vshock0,esw_frac0
            endfor
    endfor
endfor

return
end 

;-------------------------------------
;Finds the count rate map of any energy band (ebandindex) and writes
;it to a data file

;pro countrate,eband,eres,modeler,dist ; original w/o cooling depth
pro countrate,eband,eres,modeler,dist, lc, cool, vx, vy, vz, a_z, np_stream, tp_stream, vr_stream, vt_stream, ntr_frac0, nref_frac0, nen_frac0, etr_frac0, eref_frac0, een_frac0, em_frac0, kappa_ref0,vshock0,esw_frac0

;Inputs: ebandindex = 0-13, 0-7 ibex-lo, 8-13 ibex-hi
;        eres = number of energy steps (number of energies examined in
;        order to find average flux over an energy band (geometric
;        factor considered constant over one energy band))
;        modeler, the particular modeler (by name)
;        dist, the particular distribution (by name)


common distributions,dist_names,num_dist
common model_names, models, num_models, pui_model,vasyliunas,ribbon,gdf,hinterp,ntr_frac,nref_frac,ttr_frac,tref_frac,colde,esw_frac,etr_frac,eref_frac,moscow,moscow_temp,giacalone
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen
common all_sensors,Eminarr,Emaxarr,Garr,Ecenarr
common path_names,countrate_path,averagefluxplots_path,countrateplots_path,count_input_path,count_output_path
common kappa_value,kappa
common sensor_params,hi,lo,wt9,wt10,wt11,wt12,wt13,ultra

if modeler ne 'Opher_multiion_withneutrals' then countrate_map=fltarr(ntheta,nphi) else countrate_map=fltarr(ntheta,nphi,5) ; Malama MI

Energy=Eminarr[Eband]
Denergy=(Emaxarr[Eband]-Eminarr[Eband])/Eres
flux_map_test=0.

if Eband eq 9 then wt=wt9/total(wt9) else if Eband eq 10 then wt=wt10/total(wt10) else if Eband eq 11 then wt=wt11/total(wt11) else if Eband eq 12 then wt=wt12/total(wt12) else if Eband eq 13 then wt=wt13/total(wt13) 

if Eband ge 14 and Eband le 20 then Energy=Ecenarr[Eband] ; for INCA energy transmission

if modeler eq 'Opher' then opherfile
if modeler eq 'Opher_withneutrals' then opherwithneutralsfile
if modeler eq 'Opher_multiion_withneutrals' then opherwithneutralsfile ;ophermultiionwithneutralsfile
while Energy lt Emaxarr[Eband]-1e-6 do begin ; 1e-6 added by MZK since code was running energy=emaxarr[eband]
    if Eband ge 14 and Eband le 20 then begin ; for only integrating at central INCA energy
	    if pui_model eq 'Malama' then make_flux,energy,dist,modeler,flux_map1, lc, cool, vx, vy, vz, a_z, np_stream, tp_stream, vr_stream, vt_stream,eband, $
                                          ntr_frac0, nref_frac0, nen_frac0, etr_frac0, eref_frac0, een_frac0, em_frac0, kappa_ref0,vshock0,esw_frac0
            if pui_model eq 'Zirnstein' then make_flux,energy,dist,modeler,flux_map1, lc, cool, vx, vy, vz, a_z, np_stream, tp_stream, vr_stream, vt_stream,eband, $
                                          ntr_frac0, nref_frac0, nen_frac0, etr_frac0, eref_frac0, een_frac0, em_frac0, kappa_ref0,vshock0,esw_frac0
            print, 'energyband de= ', eband
            countrate_map=flux_map1*Energy*Garr[eband]
            Energy=Emaxarr[Eband] ; set energy to max energy to break while loop
    endif else begin
    	if pui_model eq 'Malama' then make_flux,energy,dist,modeler,flux_map1, lc, cool, vx, vy, vz, a_z, np_stream, tp_stream, vr_stream, vt_stream,eband, $
                                      ntr_frac0, nref_frac0, nen_frac0, etr_frac0, eref_frac0, een_frac0, em_frac0, kappa_ref0,vshock0,esw_frac0
    	if pui_model eq 'Zirnstein' then make_flux,energy,dist,modeler,flux_map1, lc, cool, vx, vy, vz, a_z, np_stream, tp_stream, vr_stream, vt_stream,eband, $
                                         ntr_frac0, nref_frac0, nen_frac0, etr_frac0, eref_frac0, een_frac0, em_frac0, kappa_ref0,vshock0,esw_frac0
    	fracfinished=(energy-Eminarr[eband])/(Emaxarr[eband]-Eminarr[eband])
    	print, energy, Eminarr[eband],Emaxarr[eband]
    	print, 'energyband de= ', eband, ' finished ',fracfinished*100.,'%'
	ect=round((energy-eminarr[eband])/denergy)
	if Eband eq 9 or Eband eq 10 or Eband eq 11 or Eband eq 12 or Eband eq 13 then countrate_map+=flux_map1*Energy*Garr[eband]*wt[ect] $
        else countrate_map+=flux_map1*Energy*Garr[eband]*denergy/(Emaxarr[Eband]-Eminarr[Eband])
    	Energy+=denergy
    endelse
endwhile

close,1
if modeler ne 'Opher_multiion_withneutrals' then begin
	openw,1,countrate_path+count_output_path+modeler+'CountRate'+'Kappa'+string(kappa,format='(f4.2)')+'Eband'+string(eband,format='(i2.2)')+'.dat'
	for i =0, ntheta-1 do begin	
		for j=0, nphi-1 do begin
	        	printf,1,countrate_map[i,j] 
    		endfor
	endfor
	close,1
endif else begin
tag = ' ' 
for k = 0,4 do begin ; Malama MI
        if k eq 0 then tag = 'PUI1' & if k eq 1 then tag = 'PUI2' & if k eq 2 then tag = 'PUI3' & if k eq 3 then tag = 'SW' & if k eq 4 then tag = 'CombPUIs'
        openw,1,countrate_path+count_output_path+modeler+tag+'CountRate'+'Kappa'+string(kappa,format='(f4.2)')+'Eband'+string(eband,format='(i2.2)')+'.dat'
       for i =0, ntheta-1 do begin     
               for j=0, nphi-1 do begin
                       printf,1,countrate_map[i,j,k] 
               endfor
       endfor
       close,1


endfor
endelse

end

;----------------------------------_
;Purpose: Reads in Opher data from secondary data files written
pro opherwithneutralsfile
;Stores in common block the  data array of proton number density,
;plasma temperature,radial velocity, tangential velocity, theta,phi,r
common secondary_path_names,secondary_path, opher_secondary_file, heerikhuisen_secondary_file, heerikhuisen_secondarynh_file, opher_withneutrals_secondary_file, opher_withneutrals_secondarynh_file, opher_twoplasmas_secondaryp1_file,opher_twoplasmas_secondaryp2_file, opher_twoplasmas_secondarynh_file, opher_withneutrals_secondarystream_file,opher_withneutrals_secondary_file1, opher_withneutrals_secondarynh_file1,opher_withneutrals_secondary_file2, opher_withneutrals_secondarynh_file2,opher_withneutrals_secondaryinterp_file,input_path,output_path
common data_opher,dataopher,sdataopher
common grid_opher,gridopher
common nhdata_opher,nhopher,snhopher ;used when neutrals
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen,r_lower,r_upper
common model_names, models, num_models, pui_model,vasyliunas,ribbon,gdf,hinterp,ntr_frac,nref_frac,ttr_frac,tref_frac,colde,esw_frac,etr_frac,eref_frac,moscow,moscow_temp, giacalone

plasmafile = secondary_path+input_path+opher_withneutrals_secondary_file
neutralfile = secondary_path+input_path+opher_withneutrals_secondarynh_file
nr=gridopher.r.num; number of r bins

i=0d
theta=0.
phi=0.
r=0.

get_lun,lun
openr,lun,plasmafile
vec=dblarr(12)
num = 0.
while not eof(lun) do begin
    readf,lun,vec
    np=vec[0]
    temp=vec[1]
    vr=vec[2]
    vtheta=vec[3]
    vphi=vec[4]
    vt=vec[5]
    i=vec[6]
    j=vec[7]
    k=vec[8]
    theta=vec[9]
    phi=vec[10]
    r=vec[11]

    sdataopher.vr[i,j,k]=vr
    sdataopher.vtheta[i,j,k]=vtheta
    sdataopher.vphi[i,j,k]=vphi
    sdataopher.vt[i,j,k]=vt
    sdataopher.density[i,j,k]=np
    if moscow_temp eq 1 then sdataopher.temp[i,j,k]=temp else sdataopher.temp[i,j,k]=temp/2.
    sdataopher.phi[i,j,k]=phi
    sdataopher.theta[i,j,k]=theta
    sdataopher.r[i,j,k]=r
    num = num +1
endwhile

close,lun
free_lun,lun
get_lun,lun
openr,lun,neutralfile

; Include Neutral Temp
vec=dblarr(8)
while not eof(lun) do begin
    readf,lun,vec
    nh=vec[0]
    th1=vec[1]
    i=vec[2]
    j=vec[3]
    k=vec[4]
    theta=vec[5]
    phi=vec[6]
    r=vec[7]
    snhopher.density[i,j,k]=nh
    snhopher.temp[i,j,k]=th1
endwhile

    if floor(sdataopher.phi[i,j,k]) ne floor(phi) then stop & if floor(sdataopher.theta[i,j,k]) ne floor(theta) then stop & if floor(sdataopher.r[i,j,k]) ne floor(r) then stop ; in previous code replaced theta,phi, r
close,lun
free_lun,lun
end

;----------------------------------------
pro make_flux,energy,dist,modeler,flux,lc, cool, vx, vy, vz, a_z, np_stream, tp_stream, vr_stream, vt_stream, eband, ntr_frac0, nref_frac0, nen_frac0, etr_frac0, eref_frac0, een_frac0, em_frac0, kappa_ref0,vshock0,esw_frac0
;Inputs:
; Energy in keV
; dist = distribution name
;Output: 
; flux: ena flux array based on merav opher's data for
; specificed energy and distribution

common secondary_path_names,secondary_path ; added by MZK for writing test boundary file
common constants_cgs, mp,kb,kev_erg, AU_cm
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen,r_lower,r_upper,r_inner,moscow_regs
common rlimits, r_limits_lower, r_limits_upper, boundaries_arr ; boundaries_arr added by MZK
common kappa_value,kappa
common magnetic_field, ophermagnetic
common radial_limits, radial_range ; ADDED BY MZK
common kappa_value,kappa,kappaPUI,kappaRef,kappaISM,kappaOH ; ADDED BY MZK
common include_survival, flag_survival, extinction ; ADDED BY MZK
common model_names, models, num_models, pui_model,vasyliunas,ribbon,gdf,hinterp,ntr_frac,nref_frac,ttr_frac,tref_frac,colde,esw_frac,etr_frac,eref_frac,moscow,moscow_temp,giacalone ; ADDED BY MZK
common sensor_params, hi,lo, wt9, wt10, wt11, wt12, wt13, ultra 

if modeler eq 'Opher' then begin 
   common grid_opher,gridopher
   common data_opher,dataopher,sdataopher	 
   nr=gridopher.r.num; number of r bins
   data= sdataopher
endif
if modeler eq 'Opher_withneutrals' or modeler eq 'Opher_multiion_withneutrals' then begin
    common grid_opher,gridopher
    common nhdata_opher,nhopher,snhopher 
    common data_opher,dataopher,sdataopher
    nr=gridopher.r.num   ; number of r bins
    data = sdataopher
    data_neutral = snhopher
endif
; Sets up arrays of neutral density (nn), plasma density (np),
; temperature (tp), theta

if modeler ne 'Opher_multiion_withneutrals' then flux=dblarr(ntheta,nphi) else flux=dblarr(ntheta,nphi,5); flux=dblarr(ntheta,nphi,4)

; initializing variables for flux calculations
nn=0D ; neutral density
np=0D ; plasma density
tp=0D ; plasma temperature
theta=0D ; heliolatitude
phi=0D ; heliolongitude
r=0D ; radius
ur=0D ; radial velocity
ut=0D ; non-radial velocity   

; converting spherical coordintes to cartesian coordinates for Tecplot output
r_xyz=dblarr(ntheta+1,nphi,nr)
theta_xyz=dblarr(ntheta+1,nphi,nr)
phi_xyz=dblarr(ntheta+1,nphi,nr)

r_xyz[0:ntheta-1,0:nphi-1,0:nr-1]=data.r
theta_xyz[0:ntheta-1,0:nphi-1,0:nr-1]=data.theta
phi_xyz[0:ntheta-1,0:nphi-1,0:nr-1]=data.phi

r_xyz[ntheta,0:nphi-1,0:nr-1]=data.r[ntheta-1,0:nphi-1,0:nr-1]
theta_xyz[ntheta,0:nphi-1,0:nr-1]=180.
phi_xyz[ntheta,0:nphi-1,0:nr-1]=data.phi[ntheta-1,0:nphi-1,0:nr-1]

x0=r_xyz*cos(phi_xyz*!pi/180.)*sin(theta_xyz*!pi/180.)
y0=r_xyz*sin(phi_xyz*!pi/180.)*sin(theta_xyz*!pi/180.)
z0=r_xyz*cos(theta_xyz*!pi/180.)

; integration limits as indices
rlower=(r_lower-r_inner)/dr
rupper=(r_upper-r_inner)/dr

; calculating energy-dependent extinction for streamlines (from read_stream.pro or read_stream_g.pro) 
vena=energy_velocity(energy)
sig_cx=xsection(energy)
ext=cool*vena*sig_cx
ext=exp(ext)
if extinction ne 1 then ext(*,*,*)=1.

; initializing relevant arrays for flux calculations
ni_rat2=make_array(ntheta,nphi,nr,/double)
np_ts=make_array(ntheta,nphi,/double)
r_ts=make_array(ntheta,nphi,/double)
r_arr=make_array(ntheta,nphi,/double)
np_rat=dblarr(ntheta,nphi)
tp_rat=dblarr(ntheta,nphi)		
dens_test1=make_array(ntheta,nphi,nr,/double)
temp_test1=make_array(ntheta,nphi,nr,/double)
dens_test2=make_array(ntheta,nphi,nr,/double)
temp_test2=make_array(ntheta,nphi,nr,/double)
dens_test3=make_array(ntheta,nphi,nr,/double)
temp_test3=make_array(ntheta,nphi,nr,/double)
dens_test_sw=make_array(ntheta,nphi,nr,/double)
temp_test_sw=make_array(ntheta,nphi,nr,/double)
temp_test_p=make_array(ntheta,nphi,nr,/double)
f_test1=make_array(ntheta,nphi,nr,/double)
lc_test=make_array(ntheta,nphi,nr)
flux_test1=make_array(ntheta,nphi,nr,/double)
flux_test2=make_array(ntheta,nphi,nr,/double)
flux_test3=make_array(ntheta,nphi,nr,/double)
flux_test_sw=make_array(ntheta,nphi,nr,/double)
rts=make_array(ntheta,nphi,/double)
rad_ts=make_array(ntheta,nphi,/double)
zirn_rat=make_array(ntheta,nphi,/double)
a=make_array(ntheta,nphi,/double)
nh_avg=make_array(ntheta,nphi,/double)
tp_ts=make_array(ntheta,nphi,/double)
v_plas=make_array(ntheta,nphi,nr,/double)
boundaries_arr0=make_array(ntheta,nphi,nr,/double)
shock_strength=make_array(ntheta,nphi,/double)
up=sqrt(data.vr^2.+data.vt^2.)

;setting up survival probability
if flag_survival eq 1 then ionization_rate,energy,bet ; reads in ionization rate
if flag_survival eq 0 then bet=dblarr(ntheta,nphi,nr)

; loop to calculate ENA flux in each grid cell and integrate along radial LOS
percentage=-1
for t=0, ntheta-1 do begin
    for p=0, nphi-1 do begin
	
	depth=where(abs(1/exp(1)-ext[t,p,*]) eq min(abs(1/exp(1)-ext[t,p,*])),dct) ; cooling depth calc
        lc[t,p]=min(depth)

	for rad=rlower,rupper do begin
		if radial_range eq 0. then begin ; FOR TS TO HP - MZK
			if boundaries_arr[t,p,rad] eq 3 or boundaries_arr[t,p,rad] eq 6 then begin; REGION 3 = HELIOSHEATH - MZK

		                if modeler eq 'Opher' then nn=.1
				if modeler eq 'Opher_withneutrals' or modeler eq 'Opher_multiion_withneutrals' then nn=data_neutral.density[t,p,rad]
             
                  ;REGULAR FLOWS USED, do not comment out inside of this portion of loop
            		       ur = data.vr[t,p,rad]*10.^5. 
                  	       ut = data.vt[t,p,rad]*10.^5.
                ;endelse 
                
               		       if ut ne ut then ut = 0D;check for NaN
               	 	       if ut gt -1.*(10.^30) THEN stupid=0 ELSE IF ut lt (10.^30) THEN stupid=0 ELSE ut=0D ;remove unphysical results
			 	theta=data.theta[t,p,rad] & phi=data.phi[t,p,rad] & r=data.r[t,p,rad]

				;SINGLE ION FLUID CASE
                 		if modeler ne 'Opher_multiion_withneutrals' then begin
              				np=data.density[t,p,rad]
                        		tp=data.temp[t,p,rad]        
                        		if boundaries_arr[t,p,rad] eq 3 then kappa0=kappa
					if boundaries_arr[t,p,rad] eq 6 then kappa0=kappaOH ; open heliosheath kappa - MZK
					if np eq 0 or tp eq 0 then stop

flux[t,p]=calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0,bet[t,p,rad])+flux[t,p] ; kappa0 and bet by MZK  

                 		endif else begin ;MULTIION FLUID CASE ------------------------------------------------------------------------

                        		np=data.density[t,p,rad]
                			tp=data.temp[t,p,rad]
   					nsw=data.density_sw[t,p,rad]
 					npt=data.density_pt[t,p,rad]
					npr=data.density_pr[t,p,rad]

                          ;TEMPERATURES FROM MODEL RESULTS USED, do not comment out inside of this portion of loop
                          		tsw  =  data.temp_sw[t,p,rad] 
					tpt = data.temp_pt[t,p,rad]
					tpr = data.temp_pr[t,p,rad]
					if boundaries_arr[t,p,rad] eq 3. then kappa0=kappa 
					if boundaries_arr[t,p,rad] eq 6. then kappa0=kappaOH

                        ;This is used for making countrate files. For regular maps make sure flows and temperatures are just as pulled from data file. Do not change.  
					flux[t,p,0]=calc_diff_flux(dist,energy,r,tsw,ur,ut,nn,nsw,theta,kappa0,bet[t,p,rad])+flux[t,p,0] ; total ion, kappa0 by MZK
                                        flux[t,p,1]=calc_diff_flux(dist,energy,r,tpt,ur,ut,nn,npt,theta,kappa0,bet[t,p,rad])+flux[t,p,1] ; SW component, kappa0 by MZK
                                        flux[t,p,2]= calc_diff_flux(dist,energy,r,tpr,ur,ut,nn,npr,theta,kappa0,bet[t,p,rad])+flux[t,p,2] ; PUI component, kappa0 by MZK

                       ;-------------------
				endelse
			endif ; ADDED BY MZK FOR REGION IF STATEMENT
		endif ; ADDED BY MZK FOR RADIAL RANGE IF STATEMENT

;\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
;/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    RADIAL RANGE 1      \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
;\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

                if radial_range eq 1. then begin ; FOR HP to 1500 AU - MZK
                        if boundaries_arr[t,p,rad] eq 2 or boundaries_arr[t,p,rad] eq 1 then begin ; REGION 2 = BEYOND HP - MZK

                                if modeler eq 'Opher' then nn=.1
                                if modeler eq 'Opher_withneutrals' or modeler eq 'Opher_multiion_withneutrals' then nn=data_neutral.density[t,p,rad]


                  ;REGULAR FLOWS USED, do not comment out inside of this portion of loop
                               ur = data.vr[t,p,rad]*10.^5.
                               ut = data.vt[t,p,rad]*10.^5.
                ;endelse 

                               if ut ne ut then ut = 0D;check for NaN
                               if ut gt -1.*(10.^30) THEN stupid=0 ELSE IF ut lt (10.^30) THEN stupid=0 ELSE ut=0D ;remove unphysical results
                                theta=data.theta[t,p,rad] & phi=data.phi[t,p,rad] & r=data.r[t,p,rad]

                                ;SINGLE ION FLUID CASE
                                if modeler ne 'Opher_multiion_withneutrals' then begin
                                        np=data.density[t,p,rad]
                                        tp=data.temp[t,p,rad]
					kappa0=kappaISM
                                        if np eq 0 or tp eq 0 then stop

flux[t,p]=calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0,bet[t,p,rad])+flux[t,p] ; kappa0 and bet by MZK  
radial_profile.ena_flux(0,t,p,rad) = calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0,bet[t,p,rad]) ; trying to find single ion flux, kappa0 and bet by MZK

                                endif else begin ;MULTIION FLUID CASE ------------------------------------------------------------------------

                                        np=data.density[t,p,rad]
                                        tp=data.temp[t,p,rad]
                                        nsw=data.density_sw[t,p,rad]/2.
                                        npu3=data.density_pu3[t,p,rad]

                          ;TEMPERATURES FROM MODEL RESULTS USED, do not comment out inside of this portion of loop
                                        tsw  =  data.temp_sw[t,p,rad]
                                        tpu3 =  data.temp_pu3[t,p,rad]
					kappa0=kappaISM
                       ;endelse 

                                endelse
			endif
		endif

;<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
;*************************************************       RADIAL RANGE 2       *************************************************
;<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

;		phi=data.phi[t,p,rad]
                if radial_range eq 2. then begin ; FOR TS to 1500 AU - MZK
			if boundaries_arr[t,p,rad] eq 3 or boundaries_arr[t,p,rad] eq 2 or boundaries_arr[t,p,rad] eq 1 or boundaries_arr[t,p,rad] eq 6 or boundaries_arr[t,p,rad] eq 4 then begin ; REGION 3 = TS TO BEYOND- MZK

                                if modeler eq 'Opher' then nn=.1
                                if modeler eq 'Opher_withneutrals' or modeler eq 'Opher_multiion_withneutrals' then nn=data_neutral.density[t,p,rad]

                  ;REGULAR FLOWS USED, do not comment out inside of this portion of loop
                               ur = data.vr[t,p,rad]*10.^5.
                               ut = data.vt[t,p,rad]*10.^5.

                               if ut ne ut then ut = 0D;check for NaN
                               if ut gt -1.*(10.^30) THEN stupid=0 ELSE IF ut lt (10.^30) THEN stupid=0 ELSE ut=0D ;remove unphysical results
                                theta=data.theta[t,p,rad] & phi=data.phi[t,p,rad] & r=data.r[t,p,rad]

                                ;SINGLE ION FLUID CASE
                                if modeler ne 'Opher_multiion_withneutrals' then begin
                                        np=data.density[t,p,rad]
                                        tp=data.temp[t,p,rad]
					if boundaries_arr[t,p,rad] eq 3. then kappa0=kappa
					if boundaries_arr[t,p,rad] eq 6. then kappa0=kappaOH; open heliosheath kappa - MZK
					if boundaries_arr[t,p,rad] lt 3. then kappa0=kappaISM; ISM kappa - MZK
					if boundaries_arr[t,p,rad] eq 4 or boundaries_arr[t,p,rad] eq 5 then kappa0=1.0
					if np eq 0 or tp eq 0 then stop

flux[t,p]=calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0,bet[t,p,rad])+flux[t,p] ; kappa0 and bet by MZK  

                                endif else begin ;MULTIION FLUID CASE ------------------------------------------------------------------------

					if boundaries_arr[t,p,rad] eq 3. then kappa0=kappa
					if boundaries_arr[t,p,rad] eq 6. then kappa0=kappaOH
                                        if boundaries_arr[t,p,rad] lt 3. then kappa0=kappaISM; Setting different kappa value for ISM - MZK
                                        if boundaries_arr[t,p,rad] eq 4 or boundaries_arr[t,p,rad] eq 5 then kappa0=1.0
					if giacalone ne 1 or boundaries_arr[t,p,rad] ne 3 then kappaEn=kappa

					if ext[t,p,rad] ne ext[t,p,rad] or ext[t,p,rad] eq 0. then ext[t,p,rad]=1.e-10

					np=data.density[t,p,rad]
					tp=data.temp[t,p,rad]
					nn=data_neutral.density[t,p,rad]
					ur=data.vr[t,p,rad]*1.e5
					ut=data.vt[t,p,rad]*1.e5

if boundaries_arr[t,p,rad] eq 3 and r le 96. then begin
	np_rat[t,p]=np/data.density[15,33,rad]
	tp_rat[t,p]=tp/data.temp[15,33,rad]
endif

nn=data_neutral.density[t,p,rad]

b=0.08
if boundaries_arr[t,p,rad] eq 4 and r gt 40 then begin 
	if p ge 15*(6./dphi) and p lt 45*(6./dphi) then np_ts[t,p]=0.2 else np_ts[t,p]=0.3
	if r lt 150 and boundaries_arr[t,p,rad-1] ne 3 and boundaries_arr[t,p,rad+1] eq 3 then rad_ts[t,p]=rad+1
        if r lt 150 and boundaries_arr[t,p,rad-1] ne 3 and boundaries_arr[t,p,rad+1] ne 3 and boundaries_arr[t,p,rad+2] eq 3 then rad_ts[t,p]=rad+2
	u1au=4.e7 ; cm/s
        n1au=8. ; cc
        nu1au=8.e-8 ; 1/s
	swe=0.5*u1au*u1au*mp/kev_erg
	sig_cx_sw=xsection(swe)
	nh_avg_ct=where(boundaries_arr[t,p,*] eq 4 and r lt 150)
	nh_avg[t,p]=mean(data_neutral.density[t,p,nh_avg_ct])
        if r lt 150 and boundaries_arr[t,p,rad-1] ne 3 then rts[t,p]=r*au_cm
	if r lt 150 and boundaries_arr[t,p,rad-1] ne 3 then a[t,p]=rts[t,p]*nh_avg[t,p]*(nu1au+sig_cx_sw*u1au*n1au)/(u1au*n1au)
endif

if boundaries_arr[t,p,rad] eq 3 and boundaries_arr[t,p,rad-1] eq 6 and a_z[t,p,rad] eq 0 and vasyliunas eq 1 then boundaries_arr[t,p,rad]=6
if boundaries_arr[t,p,rad] eq 3 and boundaries_arr[t,p,rad-1] eq 2 and a_z[t,p,rad] eq 0 and vasyliunas eq 1 then boundaries_arr[t,p,rad]=2
if boundaries_arr[t,p,rad] eq 3 and boundaries_arr[t,p,rad-1] ne 2 and boundaries_arr[t,p,rad-1] ne 1 or boundaries_arr[t,p,rad] eq 6 and boundaries_arr[t,p,rad-1] ne 2 and boundaries_arr[t,p,rad-1] ne 1 then begin

if vasyliunas eq 1 then ni_rat2[t,p,rad]=a_z[t,p,rad]/a_z[ntheta/2.,nphi/2.,ceil(r_ts[ntheta/2.,nphi/2.]-r_lower)/2] ; fixed ratios at TS that propagate along streamlines -- BEST

if giacalone eq 1 then begin
	ntr_frac=ntr_frac0[t,p,rad]
	nref_frac=nref_frac0[t,p,rad]
	nen_frac=nen_frac0[t,p,rad]
	etr_frac=etr_frac0[t,p,rad]
	eref_frac=eref_frac0[t,p,rad]
	een_frac=een_frac0[t,p,rad]
        em_frac=em_frac0[t,p,rad]
	kappaEn=kappa_ref0[t,p,rad]
        vshock=vshock0[t,p,rad]
	esw_frac=esw_frac0[t,p,rad]

        nsw_frac=1-ntr_frac-nref_frac

        ttr_frac=etr_frac/ntr_frac
        tref_frac=eref_frac/(nref_frac*(1-nen_frac))
        tsw_frac=0.04/nsw_frac
        ten_frac=(1-etr_frac-eref_frac-0.04)/(nen_frac*nref_frac)
        if colde eq 1 and giacalone eq 1 then tp=2*tp/(1+tsw_frac)
endif

if finite(ni_rat2[t,p,rad]) ne 1 then stop
zirn_rat[t,p]=a[t,p]/a[ntheta/2.,nphi/2.]
                ; Primary Method for Modeling PUIs
                if vasyliunas eq 1 then np1=ni_rat2[t,p,rad]*ntr_frac*np*ext[t,p,rad] else np1=ntr_frac*np*ext[t,p,rad] ; MHD Profile w/ Vasyliunas scaling
		if vasyliunas eq 1 then np2=ni_rat2[t,p,rad]*nref_frac*np*ext[t,p,rad] else np2=nref_frac*(1-nen_frac)*np*ext[t,p,rad] ; MHD Profile w/ Vasyliunas scaling
		if giacalone ne 1 then nsw=(1-ni_rat2[t,p,rad]*(ntr_frac+nref_frac))*np*ext[t,p,rad] else nsw=nsw_frac*np*ext[t,p,rad] ; solar wind
                if giacalone ne 1 then np3=np-nsw-np1-np2 else np3=nref_frac*nen_frac*np*ext[t,p,rad]; np3=nen_frac*np*ext[t,p,rad]; Locally created PUIs
                if np3 le 0 then np3=1e-10
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then np3=np3 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then np3=1e-10
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then np1=np1*0.8 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then np1=1e-10
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then np2=np2 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then np2=1e-10
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then nsw=np1*0.2 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then nsw=1e-10
		
                if boundaries_arr[t,p,rad] eq 6 then ur=ur;30e5
		if boundaries_arr[t,p,rad] eq 6 then ut=ut;0.

                if colde eq 1 and giacalone ne 1 then tp=2*tp/(1+((1-ttr_frac*ntr_frac-tref_frac*nref_frac)/(1-ni_rat2[t,p,rad]*(ntr_frac+nref_frac))))
		if giacalone ne 1 then tp1=ttr_frac*tp/ni_rat2[t,p,rad] else tp1=tp*ttr_frac ; Transmitted PUIs from interstellar neutrals
                if giacalone ne 1 then tp2=tref_frac*tp/ni_rat2[t,p,rad] else tp2=tp*tref_frac ; Reflected PUIs
                if giacalone ne 1 then tp3=1e-10 else tp3=tp*ten_frac
		if giacalone ne 1 then tsw=tp*(1-ttr_frac*ntr_frac-tref_frac*nref_frac)/(1-ni_rat2[t,p,rad]*(ntr_frac+nref_frac)) else tsw=tp*tsw_frac
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tp3=tp3 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tp3=1e-10;3.*tp
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tp1=tp1 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tp1=1e-10; 1e-10
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tp2=tp2 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tp2=1e-10; 1e-10
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tsw=(1+(bsw/0.15)^2.)*tp1 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tsw=1e-10; 1e-10

if tp1 le 0 or np1 le 0 or tp2 le 0 or np2 le 0 or tp3 le 0 or np3 le 0 or tsw le 0 or nsw le 0 then stop
if tp1 ne tp1 or np1 ne np1 or tp2 ne tp2 or np2 ne np2 or tp3 ne tp3 or np3 ne np3 or tsw ne tsw or nsw ne nsw then stop

dens_test1[t,p,rad]=np1
temp_test1[t,p,rad]=tp1
dens_test2[t,p,rad]=np2
temp_test2[t,p,rad]=tp2
dens_test3[t,p,rad]=np3
temp_test3[t,p,rad]=tp3
dens_test_sw[t,p,rad]=nsw
temp_test_sw[t,p,rad]=tsw
temp_test_p[t,p,rad]=tp

data.vr[t,p,rad]=ur
data.vt[t,p,rad]=ut

v_plasma = plasma_frame(energy,ur,ut)
f_test1[t,p,rad]=double(maxwellian_fn(tp1,np1,v_plasma))
v_plas[t,p,rad]=v_plasma
if rad eq lc[t,p] then lc_test[t,p,rad]=1. else lc_test[t,p,rad]=0
if boundaries_arr[t,p,rad] eq 3 and boundaries_arr[t,p,rad-1] eq 4 then shock_strength[t,p]=up[t,p,rad-2]/up[t,p,rad]

if ur ne ur or ut ne ut or energy ne energy then stop
                                   if boundaries_arr[t,p,rad] eq 3 then begin
					flux[t,p,0]=calc_diff_flux(dist,energy,r,tp1,ur,ut,nn,np1,theta,kappaPUI,bet[t,p,rad])+flux[t,p,0] ; Trans PUI, kappa0 by MZK
					flux[t,p,1]=calc_diff_flux(dist,energy,r,tp2,ur,ut,nn,np2,theta,kappaRef,bet[t,p,rad])+flux[t,p,1] ; Refl PUI (ENA), kappa0 by MZK
					if giacalone ne 1 then flux[t,p,2]=calc_diff_flux(dist,energy,r,tp3,ur,ut,nn,np3,theta,kappaEn,bet[t,p,rad])+flux[t,p,2] else $
                                                               flux[t,p,2]=calc_diff_flux_zank(dist,energy,r,ur,ut,nn,np3,theta,kappaEn,bet[t,p,rad],een_frac,em_frac,vshock)+flux[t,p,2] ; local PUI, kappa0 by MZK
					flux[t,p,3]=calc_diff_flux(dist,energy,r,tsw,ur,ut,nn,nsw,theta,kappa0,bet[t,p,rad])+flux[t,p,3] ; SW ion, kappa0 by MZK
                                    endif else if boundaries_arr[t,p,rad] eq 6 then begin
                                        flux[t,p,0]=calc_diff_flux(dist,energy,r,tp1,ur,ut,nn,np1,theta,kappaPUI,bet[t,p,rad])+flux[t,p,0] ; Trans PUI, kappa0 by MZK
                                        flux[t,p,1]=calc_diff_flux(dist,energy,r,tp2,ur,ut,nn,np2,theta,kappaRef,bet[t,p,rad])+flux[t,p,1] ; Refl PUI (ENA), kappa0 by MZK
                                        if giacalone ne 1 then flux[t,p,2]=calc_diff_flux(dist,energy,r,tp3,ur,ut,nn,np3,theta,kappaEn,bet[t,p,rad])+flux[t,p,2] else $
                                                               flux[t,p,2]=calc_diff_flux_zank(dist,energy,r,ur,ut,nn,np3,theta,kappaEn,bet[t,p,rad],een_frac,em_frac,vshock)+flux[t,p,2] ; local PUI, kappa0 by MZK
                                        flux[t,p,3]=calc_diff_flux(dist,energy,r,tsw,ur,ut,nn,nsw,theta,kappa0,bet[t,p,rad])+flux[t,p,3] ; SW ion, kappa0 by MZK
                                    endif
 
					flux_test1[t,p,rad]=calc_diff_flux(dist,energy,r,tp1,ur,ut,nn,np1,theta,kappaPUI,bet[t,p,rad])
					if boundaries_arr[t,p,rad] eq 6 then flux_test2[t,p,rad]=calc_diff_flux(dist,energy,r,tp2,ur,ut,nn,np2,theta,kappaOH,bet[t,p,rad])
					if boundaries_arr[t,p,rad] ne 6 then flux_test2[t,p,rad]=calc_diff_flux(dist,energy,r,tp2,ur,ut,nn,np2,theta,kappaRef,bet[t,p,rad])
					if giacalone ne 1 then flux_test3[t,p,rad]=calc_diff_flux(dist,energy,r,tp3,ur,ut,nn,np3,theta,kappaEn,bet[t,p,rad]) else $
                                                               flux_test3[t,p,rad]=calc_diff_flux_zank(dist,energy,r,ur,ut,nn,np3,theta,kappaEn,bet[t,p,rad],een_frac,em_frac,vshock)
					flux_test_sw[t,p,rad]=calc_diff_flux(dist,energy,r,tsw,ur,ut,nn,nsw,theta,kappa0,bet[t,p,rad])

endif
                                endelse
                        endif
                endif

;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   	endfor
if modeler eq 'Opher_multiion_withneutrals' then flux[t,p,4]=flux[t,p,0]+flux[t,p,1]+flux[t,p,2]+flux[t,p,3]
    endfor	
endfor

if giacalone eq 1 then begin ; repopulating some arrays for plotting purposes
a_z=nen_frac0 ; accelerated fraction
ni_rat2=een_frac0 ; min acc. energy
np_stream=em_frac0 ; max acc. energy
tp_stream=kappa_ref0 ; delta
vr_stream=vshock0 ; speed of shock
endif

close,1

if energy gt 1.1 and energy lt 1.14 then begin
   openw,1,'GH_1_11keV.dat'
   for i=30,60 do printf, 1, data.r[19,29,i],dens_test1[19,29,i],temp_test1[19,29,i],flux_test1[19,29,i],dens_test2[19,29,i],temp_test2[19,29,i],flux_test2[19,29,i],dens_test3[19,29,i],temp_test3[19,29,i],flux_test3[19,29,i] 
   close,1
endif

if energy gt 4.1 and energy lt 4.5 then begin
   openw,1,'GH_4_29keV.dat'
   for i=30,60 do printf, 1, data.r[19,29,i],dens_test1[19,29,i],temp_test1[19,29,i],flux_test1[19,29,i],dens_test2[19,29,i],temp_test2[19,29,i],flux_test2[19,29,i],dens_test3[19,29,i],temp_test3[19,29,i],flux_test3[19,29,i]
   close,1
endif

emin=[0.68,1.1,1.7,2.65,4.1,8.1,17.0,28.5,42.8,33.,55.,102.,168.]
emax=[0.73,1.14,1.8,2.75,4.5,8.7,18.5,30.0,44.8,34.,57.,106.,174.]

; Making output files for each energy band
if energy gt emin[eband-9] and energy lt emax[eband-9] then begin
print, "GENERATING OUTPUT FILE FOR EBAND ", eband ," AT ENERGY ", energy 
openw,1,secondary_path+'output/tec/denstemptest_'+'Eband'+string(eband,format='(i2.2)')+'.dat'
for i=0,ntheta-1 do begin
        for j=0, nphi-1 do begin
                for k=rlower, rupper-1 do begin
                        if ext[i,j,k] ne ext[i,j,k] then ext[i,j,k]=0.
                        printf,1,x0[i,j,k],y0[i,j,k],z0[i,j,k],boundaries_arr[i,j,k],i,j,k,data.theta[i,j,k],data.phi[i,j,k],data.r[i,j,k], dens_test1[i,j,k],temp_test1[i,j,k], dens_test2[i,j,k],temp_test2[i,j,k], dens_test3[i,j,k],temp_test3[i,j,k], dens_test_sw[i,j,k],temp_test_sw[i,j,k],data.density[i,j,k],temp_test_p[i,j,k],data.vr[i,j,k], data.vt[i,j,k],f_test1[i,j,k],lc_test[i,j,k],ext[i,j,k], vx[i,j,k], vy[i,j,k], vz[i,j,k],flux_test1[i,j,k],flux_test2[i,j,k],flux_test3[i,j,k],flux_test_sw[i,j,k], ophermagnetic.bx[i,j,k], ophermagnetic.by[i,j,k], ophermagnetic.bz[i,j,k], data_neutral.density[i,j,k], v_plas[i,j,k],a_z[i,j,k],ni_rat2[i,j,k],np_stream[i,j,k], tp_stream[i,j,k], vr_stream[i,j,k], vt_stream[i,j,k]
                endfor
        endfor
endfor
i=ntheta-1
for j=0,nphi-1 do begin
        for k=rlower,rupper-1 do begin
                if ext[i,j,k] ne ext[i,j,k] then ext[i,j,k]=0.
                        printf,1,x0[i+1,j,k],y0[i+1,j,k],z0[i+1,j,k],boundaries_arr[i,j,k],i,j,k,theta_xyz[i+1,j,k],data.phi[i,j,k],data.r[i,j,k], dens_test1[i,j,k],temp_test1[i,j,k], dens_test2[i,j,k],temp_test2[i,j,k], dens_test3[i,j,k],temp_test3[i,j,k], dens_test_sw[i,j,k],temp_test_sw[i,j,k],data.density[i,j,k],temp_test_p[i,j,k],data.vr[i,j,k], data.vt[i,j,k],f_test1[i,j,k],lc_test[i,j,k],ext[i,j,k], vx[i,j,k], vy[i,j,k], vz[i,j,k],flux_test1[i,j,k],flux_test2[i,j,k],flux_test3[i,j,k],flux_test_sw[i,j,k], ophermagnetic.bx[i,j,k], ophermagnetic.by[i,j,k], ophermagnetic.bz[i,j,k], data_neutral.density[i,j,k], v_plas[i,j,k],a_z[i,j,k],ni_rat2[i,j,k],np_stream[i,j,k], tp_stream[i,j,k], vr_stream[i,j,k], vt_stream[i,j,k]
        endfor
endfor
close,1
print, "Cooling length for tail: ", data.r[22,0,lc[22,0]]
endif


return
end

;---------------------------------
function plasma_frame,E,vr,vn

;function plasma_frame, E, Tp, vr, vn

;take observed energy and calculate velocity in frame of heliosheath plasma
;vr and vn are radial and normal components of downstream solar wind flow 
common constants_cgs, mp,kb,kev_erg,AU_cm

VenaErf=SQRT(2.0*E*kev_erg/mp)   ; ion velocity in energy band (cm/s from KE = 1/2mv^2)

vplasma=SQRT((VenaErf+vr)^2.+vn^2.)

return, vplasma
end
;------------------------------
function plasma_frame2,E,vr,vn,tp
;take observed energy and calculate velocity in frame of heliosheath plasma
;vr and vn are radial and normal components of downstream solar wind flow 
common constants_cgs, mp,kb,kev_erg,AU_cm

;assume all of ENA velocity is radial
VenaErf=SQRT(2.0*E*kev_erg/mp)   ; in m/s from KE = 1/2mv^2
vplasma=SQRT((VenaErf+vr)^2.+vn^2.)

vth = sqrt(2.*kb*tp/mp) ; thermal velocity

w=vplasma/vth

v1=exp(-w^2.)/sqrt(!pi)
v2=(w+1./(2.*w))*erf(w)
vplasma2=vth*(v1+v2) ; relative velocity

return,vplasma2
end
;------------------------------
function ENA_frame,E,vr 
common constants_cgs, mp,kb,kev_erg,AU_cm

;assume all of ENA velocity is radial
VenaErf=SQRT(2.0*E*kev_erg/mp)   ; in m/s from KE = 1/2mv^2
vplasma=SQRT(VenaErf^2.+vr^2.)

return,vplasma
end
;------------------------------
function energy_velocity,E
;converts an energy to a velocity

common constants_cgs, mp,kb,kev_erg,AU_cm
v =SQRT(2.0*E*kev_erg/mp)
return,v

end
;-------------------------------

function vrel,T,vplasma,vr_h ; by MZK
;calculate the relative velocity between parent proton and neutral from Heerikhuisen et al. (2006)
common constants_cgs, mp,kb,kev_erg,AU_cm

vth = sqrt(2.*kb*T/mp) ; thermal velocity

w=sqrt((vplasma-vr_h)^2.)/vth 

v1=exp(-w^2.)/sqrt(!pi)
v2=(w+1./(2.*w))*erf(w)
v=vth*(v1+v2) ; relative velocity

return, v

end
;-------------------------------

pro ionization_rate,Ea,bet ; by MZK
; ionization rate integration along line-of-sight
common grid_opher,gridopher
common data_opher,dataopher,sdataopher
common nhdata_opher,nhopher,snhopher
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen,r_lower,r_upper,r_inner,moscow_regs 
common constants_cgs, mp,kb,kev_erg,AU_cm
common rlimits, r_limits_lower, r_limits_upper, boundaries_arr
nr=gridopher.r.num

bet=make_array(ntheta,nphi,nr)
vplasma=make_array(ntheta,nphi,nr)
sigma=make_array(ntheta,nphi,nr)
;sigma=xsection(Ea)
vena=energy_velocity(Ea)
E_rel=make_array(ntheta,nphi,nr)

tp=sdataopher.temp
vr=sdataopher.vr*1.e5
vn=sdataopher.vt*1.e5
np=sdataopher.density
nh=snhopher.density
beta0=6.0e-7
r1=AU_cm

for i=0,ntheta-1 do begin
	for j=0,nphi-1 do begin
		for k=(100-r_inner)/dr,nr-1 do begin ; to do survival probability into 100 AU
			if boundaries_arr[i,j,k] eq 3 then begin
				vplasma[i,j,k]=plasma_frame2(Ea,vr[i,j,k],vn[i,j,k],tp[i,j,k])
				E_rel[i,j,k]=0.5*mp*(vplasma[i,j,k])^2.
				if E_rel[i,j,k] ne E_rel[i,j,k] then stop
				sigma[i,j,k]=xsection(E_rel[i,j,k])
				bet[i,j,k]=double(sigma[i,j,k]*vplasma[i,j,k]*np[i,j,k]*dr)
				bet[i,j,k]=bet[i,j,k]+bet[i,j,k-1]	
			endif
		endfor
	endfor
endfor

end
