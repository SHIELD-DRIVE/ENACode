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

;--------------------------------------
function calc_diff_flux_multiion,dist,Ea,r,Tp,Ur,Un,nn,np,theta,Uphi,Utheta,kappa0,bet ; kappa0 by MZK
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
;Current of ENAs
common kappa_value,kappa,kappaPUI,kappaRef,kappaISM,kappaOH
common constants_cgs, mp,kb,kev_erg, AU_cm
common distributions,dist_names
;put energy into solar wind frame and convert to velocity

if kappa0 le 1.50000 and dist eq 'Kappa'+string(kappa,format='(f4.2)') then dist='Maxwellian'
if kappa0 gt 1.50000 and dist eq 'Kappa'+string(kappa,format='(f4.2)') then dist='Kappa'+string(kappa,format='(f4.2)')
v_plasma = plasma_frame(Ea,Ur,Un)
;v_plasma = plasma_frame2(Ea,Ur,Uphi,Utheta)


if dist eq 'Kappa'+string(kappa,format='(f4.2)') then f=double(kappa_fn(Tp,np,v_plasma,kappa0))
if dist eq 'Maxwellian' then f=double(maxwellian_fn(Tp,np,v_plasma))
fluxconstant=double(fluxconstants(Ea,r,nn,theta,bet)) ;Sigma*Nn*S*survival probability
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

; PARKER TRANSPORT TESTING
;parker_transport, lc, cool, vx, vy, vz, ntr_frac0, nref_frac0, nen_frac0, etr_frac0, eref_frac0, een_frac0,em_frac0, kappa_ref0,vshock0,esw_frac0
   
;if extinction ne 0 then read_stream, lc, cool, vx, vy, vz, a_z, np_stream, tp_stream, vr_stream, vt_stream
if extinction eq 1 and giacalone eq 0 or vasyliunas eq 1 then read_stream, lc, cool, vx, vy, vz, a_z, np_stream, tp_stream, vr_stream, vt_stream ; GIACALONE
;read_stream, lc, cool, vx, vy, vz, a_z, np_stream, tp_stream, vr_stream, vt_stream ; GIACALONE
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

;if modeler ne 'Opher_multiion_withneutrals' then countrate_map=fltarr(ntheta,nphi) else countrate_map=fltarr(ntheta,nphi,3) ; original MI
;if modeler ne 'Opher_multiion_withneutrals' then countrate_map=fltarr(ntheta,nphi) else countrate_map=fltarr(ntheta,nphi,4) ; Zank MI
if modeler ne 'Opher_multiion_withneutrals' then countrate_map=fltarr(ntheta,nphi) else countrate_map=fltarr(ntheta,nphi,5) ; Malama MI

Energy=Eminarr[Eband]
Denergy=(Emaxarr[Eband]-Eminarr[Eband])/Eres
flux_map_test=0.

if Eband eq 9 then wt=wt9/total(wt9) else if Eband eq 10 then wt=wt10/total(wt10) else if Eband eq 11 then wt=wt11/total(wt11) else if Eband eq 12 then wt=wt12/total(wt12) else if Eband eq 13 then wt=wt13/total(wt13) 

if Eband ge 14 and Eband le 20 then Energy=Ecenarr[Eband] ; for INCA energy transmission

if modeler eq 'Opher' then opherfile
if modeler eq 'Heerikhuisen' then heerikhuisenfile
if modeler eq 'Opher_withneutrals' then opherwithneutralsfile
if modeler eq 'Opher_multiion_withneutrals' then opherwithneutralsfile ;ophermultiionwithneutralsfile
while Energy lt Emaxarr[Eband]-1e-6 do begin ; 1e-6 added by MZK since code was running energy=emaxarr[eband]
;    make_flux,energy,dist,modeler,flux_map1 ; original w/o cooling depth
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
;See Lab book pg. 48 for equation
;    print, 'band flux at (15,30) = ', flux_map1(15,30)
	ect=round((energy-eminarr[eband])/denergy)
;	countrate_map+=flux_map1*Energy*Garr[eband]*denergy/(Emaxarr[Eband]-Eminarr[Eband])
	if Eband eq 9 or Eband eq 10 or Eband eq 11 or Eband eq 12 or Eband eq 13 then countrate_map+=flux_map1*Energy*Garr[eband]*wt[ect] $
        else countrate_map+=flux_map1*Energy*Garr[eband]*denergy/(Emaxarr[Eband]-Eminarr[Eband])
;	countrate_map_=flux_map1
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
;for k = 0,2 do begin ; original MI
;        if k eq 0 then tag = 'totplasma' & if k eq 1 then tag = 'swplasma' & if k eq 2 then tag = 'pu3plasma'
;        openw,1,countrate_path+modeler+tag+'CountRate'+'Kappa'+string(kappa,format='(f4.2)')+'Eband'+string(eband,format='(i2.2)')+'.dat'
;	for i =0, ntheta-1 do begin	
;		for j=0, nphi-1 do begin
;	        	printf,1,countrate_map[i,j,k] 
;    		endfor
;	endfor
;	close,1
;for k = 0,3 do begin ; Zank MI
;        if k eq 0 then tag = 'swplasma' & if k eq 1 then tag = 'ptplasma' & if k eq 2 then tag = 'prplasma' & if k eq 3 then tag = 'combplasma'
;        openw,1,countrate_path+modeler+tag+'CountRate'+'Kappa'+string(kappa,format='(f4.2)')+'Eband'+string(eband,format='(i2.2)')+'.dat'
;       for i =0, ntheta-1 do begin     
;               for j=0, nphi-1 do begin
;                       printf,1,countrate_map[i,j,k] 
;               endfor
;       endfor
;       close,1
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

;-----------------------------------------
;Purpose: Reads in Opher data from secondary data files written
pro opherfile
;Stores in common block the  data array of proton number density,
;plasma temperature,radial velocity, tangential velocity, theta,phi,r
common secondary_path_names,secondary_path, opher_secondary_file, heerikhuisen_secondary_file, heerikhuisen_secondarynh_file, opher_withneutrals_secondary_file, opher_withneutrals_secondarynh_file, opher_twoplasmas_secondaryp1_file,opher_twoplasmas_secondaryp2_file, opher_twoplasmas_secondarynh_file, opher_withneutrals_secondarystream_file,opher_withneutrals_secondary_file1, opher_withneutrals_secondarynh_file1,opher_withneutrals_secondary_file2, opher_withneutrals_secondarynh_file2,opher_withneutrals_secondaryinterp_file,input_path,output_path
common data_opher,dataopher,sdataopher
common grid_opher,gridopher

plasmafile = secondary_path+input_path+opher_secondary_file
nr=gridopher.r.num; number of r bins

i=0d
theta=0.
phi=0.
r=0.

close,1
openr,1,plasmafile
vec=dblarr(12)
while not eof(1) do begin
    readf,1,vec
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
    sdataopher.temp[i,j,k]=temp
    sdataopher.phi[i,j,k]=phi
    sdataopher.theta[i,j,k]=theta
    sdataopher.r[i,j,k]=r
    sdataopher.r[i,j,k]=r

endwhile

close,1

end

;----------------------------------------
;Purpose: Reads in Heerikhuisen data from secondary data files written
pro heerikhuisenfile
common secondary_path_names,secondary_path, opher_secondary_file, heerikhuisen_secondary_file, heerikhuisen_secondarynh_file, opher_withneutrals_secondary_file, opher_withneutrals_secondarynh_file, opher_twoplasmas_secondaryp1_file,opher_twoplasmas_secondaryp2_file, opher_twoplasmas_secondarynh_file, opher_withneutrals_secondarystream_file,opher_withneutrals_secondary_file1, opher_withneutrals_secondarynh_file1,opher_withneutrals_secondary_file2, opher_withneutrals_secondarynh_file2,opher_withneutrals_secondaryinterp_file,input_path,output_path
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen
common nhdata_heerikhuisen,nhheer,snhheer
common data_heerikhuisen,dataheer,sdataheer
common grid_heerikhuisen,gridheer

plasmafile = secondary_path+input_path+heerikhuisen_secondary_file
nhdfile = secondary_path+input_path+heerikhuisen_secondarynh_file

nr=gridheer.r.num; number of r bins

i=0d
theta=0.
phi=0.
r=0.

close,1
openr,1,plasmafile
vec=dblarr(12)
while not eof(1) do begin
    readf,1,vec
    np=vec[0]
    tp=vec[1]
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

    sdataheer.vr[i,j,k]=vr
    sdataheer.vtheta[i,j,k]=vtheta
    sdataheer.vphi[i,j,k]=vphi
    sdataheer.vt[i,j,k]=vt
    sdataheer.density[i,j,k]=np
    sdataheer.temp[i,j,k]=tp
    sdataheer.phi[i,j,k]=phi
    sdataheer.theta[i,j,k]=theta
    sdataheer.r[i,j,k]=r
endwhile

close,1

close,2
openr,2,nhdfile
vec2=dblarr(7)
while not eof(2) do begin
    readf,2,vec2
    nn=vec2[0]
    i=vec2[1]
    j=vec2[2]
    k=vec2[3]
    theta=vec2[4]
    phi=vec2[5]
    r=vec2[6]

    snhheer.density[i,j,k]=nn
    if sdataheer.phi[i,j,k] ne phi then stop & if sdataheer.theta[i,j,k] ne theta then stop & if sdataheer.r[i,j,k] ne r then stop ; in previous code replaced theta,phi, r
endwhile
close,2

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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Routine to rotate phi using rotation matrix from secondary_multiion_opher_withneutrals.pro.
; Can instead use phi'=phi+rot_angle, where rot_angle is some rotation angle.
; While this did not have an effect, if want to try and make work, be sure to 
; change phi, theta, and r in the definitions of sdataopher to accomodate for rotated
; angles.
;
;theta_dtor=theta*!dtor
;phi_dtor=phi*!dtor
;x_rot=r*sin(theta_dtor)*cos(phi_dtor)
;y_rot=r*sin(phi_dtor)*sin(theta_dtor)
;z_rot=r*cos(theta_dtor)
;vec_rot=[x_rot,y_rot,z_rot]
;vec_null=[0.,0.,0.]
;gamma=(255.-180.)*!dtor
;matrix=[[cos(gamma), -1.*sin(gamma), 0.],$
;        [sin(gamma),     cos(gamma), 0.],$
;        [        0.,             0., 1.]]   
;vec_temp=matrix##vec_rot
;vec_null[0]=vec_temp[0,0] & vec_null[1]=vec_temp[0,1] & vec_null[2]=vec_temp[0,2]
;r=sqrt(vec_null[0]^2.+vec_null[1]^2.+vec_null[2]^2.)
;theta=acos(vec_null[2]/r)/!dtor
;phi=atan(vec_null[1],vec_null[2])/!dtor
;if phi lt 0 then phi=360+phi ; plus since phi is negative
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


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

;openw,1, secondary_path+'inputtest.dat'
;for i=0, ntheta-1 do begin
;	for j=0, nphi-1 do begin
;		for k=0,nr-1 do begin
;			printf, 1, i, j, k, sdataopher.theta[i,j,k], sdataopher.phi[i,j,k], sdataopher.r[i,j,k], sdataopher.density[i,j,k], sdataopher.temp[i,j,k], sdataopher.vr[i,j,k], sdataopher.vt[i,j,k]
;		endfor
;	endfor
;endfor
;close,1

close,lun
free_lun,lun
get_lun,lun
openr,lun,neutralfile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Differentiating between input neutral densities
;vec=dblarr(7)
;while not eof(lun) do begin
;    readf,lun,vec
;    nh=vec[0]
;    i=vec[1]
;    j=vec[2]
;    k=vec[3]
;    theta=vec[4]
;    phi=vec[5]
;    r=vec[6]
;    snhopher.density[i,j,k]=nh
;endwhile

; Include Pop I Neutral Temp
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

; Include Pop I Neutral Temp & Weighted Neutral Speed
;vec=dblarr(10)
;while not eof(lun) do begin
;    readf,lun,vec
;    nh=vec[0]
;    th1=vec[1]
;    vr_nh=vec[2]
;    vt_nh=vec[3]
;    i=vec[4]
;    j=vec[5]
;    k=vec[6]
;    theta=vec[7]
;    phi=vec[8]
;    r=vec[9]
;    snhopher.density[i,j,k]=nh
;    snhopher.temp[i,j,k]=th1
;    snhoper.vr[i,j,k]=vr_nh
;    snhoper.vt[i,j,k]=vt_nh
;endwhile


; Include Neutral Densities and Speeds
;vec=dblarr(18)
;while not eof(lun) do begin
;readf,lun,vec
;	nh1=vec[0]
;	vr_1=vec[1]
;	vt_1=vec[2]
;	nh2=vec[3]
;	vr_2=vec[4]
;	vt_2=vec[5]
;	nh3=vec[6]
;	vr_3=vec[7]
;	vt_3=vec[8]
;	nh4=vec[9]
;	vr_4=vec[10]
;	vt_4=vec[11]
;	i=vec[12]
;	j=vec[13]
;	k=vec[14]
;	theta=vec[15]
;	phi=vec[16]
;	r=vec[17]
;	snhopher.density[i,j,k,0]=nh1
;	snhopher.density[i,j,k,1]=nh2
;	snhopher.density[i,j,k,2]=nh3
;	snhopher.density[i,j,k,3]=nh4
;	snhopher.vr[i,j,k,0]=vr_1
;	snhopher.vr[i,j,k,1]=vr_2
;	snhopher.vr[i,j,k,2]=vr_3
;	snhopher.vr[i,j,k,3]=vr_4
;	snhopher.vt[i,j,k,0]=vt_1
;	snhopher.vt[i,j,k,1]=vt_2
;	snhopher.vt[i,j,k,2]=vt_3
;	snhopher.vt[i,j,k,3]=vt_4
;endwhile
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    if floor(sdataopher.phi[i,j,k]) ne floor(phi) then stop & if floor(sdataopher.theta[i,j,k]) ne floor(theta) then stop & if floor(sdataopher.r[i,j,k]) ne floor(r) then stop ; in previous code replaced theta,phi, r
close,lun
free_lun,lun
end

;----------------------------------------
;Purpose: Reads in Opher data from secondary data files written
pro ophermultiionwithneutralsfile
;Stores in common block the  data array of proton number density,
;plasma temperature,radial velocity, tangential velocity, theta,phi,r
common secondary_path_names,secondary_path, opher_secondary_file, heerikhuisen_secondary_file, heerikhuisen_secondarynh_file, opher_withneutrals_secondary_file, opher_withneutrals_secondarynh_file, opher_twoplasmas_secondaryp1_file,opher_twoplasmas_secondaryp2_file, opher_twoplasmas_secondarynh_file, opher_withneutrals_secondarystream_file,opher_withneutrals_secondary_file1, opher_withneutrals_secondarynh_file1,opher_withneutrals_secondary_file2, opher_withneutrals_secondarynh_file2,opher_withneutrals_secondaryinterp_file,input_path,output_path
common data_opher,dataopher,sdataopher
common grid_opher,gridopher
common nhdata_opher,nhopher,snhopher ;used when neutrals

plasmafile = secondary_path+input_path+opher_withneutrals_secondary_file
neutralfile = secondary_path+input+path+opher_withneutrals_secondarynh_file
nr=gridopher.r.num; number of r bins

i=0d
theta=0.
phi=0.
r=0.

get_lun,lun
openr,lun,plasmafile
;vec=dblarr(16)
vec=dblarr(18)
num = 0.
while not eof(lun) do begin
    readf,lun,vec
;    np=vec[0]
;    temp=vec[1]
;    nsw = vec[2]
;    tempsw = vec[3]
;    npu3 = vec[4]
;    temppu3 = vec[5]
;    vr=vec[6]
;    vtheta=vec[7]
;    vphi=vec[8]
;    vt=vec[9]
;    i=vec[10]
;    j=vec[11]
;    k=vec[12]
;    theta=vec[13]
;    phi=vec[14]
;    r=vec[15]

    np=vec[0]
    temp=vec[1]
    nsw = vec[2]
    tempsw = vec[3]
    npt = vec[4]
    temppt = vec[5]
    npr = vec[6]
    temppr = vec[7]
    vr=vec[8]
    vtheta=vec[9]
    vphi=vec[10]
    vt=vec[11]
    i=vec[12]
    j=vec[13]
    k=vec[14]
    theta=vec[15]
    phi=vec[16]
    r=vec[17]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Routine to rotate phi using rotation matrix from secondary_multiion_opher_withneutrals.pro.
; Can instead use phi'=phi+rot_angle, where rot_angle is some rotation angle.
; While this did not have an effect, if want to try and make work, be sure to 
; change phi, theta, and r in the definitions of sdataopher to accomodate for rotated
; angles.
;
;theta_dtor=theta*!dtor
;phi_dtor=phi*!dtor
;x_rot=r*sin(theta_dtor)*cos(phi_dtor)
;y_rot=r*sin(phi_dtor)*sin(theta_dtor)
;z_rot=r*cos(theta_dtor)
;vec_rot=[x_rot,y_rot,z_rot]
;vec_null=[0.,0.,0.]
;gamma=(255.-180.)*!dtor
;matrix=[[cos(gamma), -1.*sin(gamma), 0.],$
;        [sin(gamma),     cos(gamma), 0.],$
;        [        0.,             0., 1.]]   
;vec_temp=matrix##vec_rot
;vec_null[0]=vec_temp[0,0] & vec_null[1]=vec_temp[0,1] & vec_null[2]=vec_temp[0,2]
;r_rot=sqrt(vec_null[0]^2.+vec_null[1]^2.+vec_null[2]^2.)
;theta=acos(vec_null[2]/r_rot)
;phi=atan(vec_null[1],vec_null[2])
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    sdataopher.vr[i,j,k]=vr
    sdataopher.vt[i,j,k]=vt
    sdataopher.density[i,j,k]=np
    sdataopher.temp[i,j,k]=temp
    sdataopher.density_sw[i,j,k]=nsw
    sdataopher.temp_sw[i,j,k]=tempsw
    sdataopher.density_pt[i,j,k]=npt
    sdataopher.temp_pt[i,j,k]=temppt
    sdataopher.density_pr[i,j,k]=npr
    sdataopher.temp_pr[i,j,k]=temppr
 ;   sdataopher.density_pu3[i,j,k]=npu3
 ;   sdataopher.temp_pu3[i,j,k]=temppu3
    sdataopher.vtheta[i,j,k] = vtheta
    sdataopher.vphi[i,j,k]= vphi
    sdataopher.phi[i,j,k]=phi
    sdataopher.theta[i,j,k]=theta
    sdataopher.r[i,j,k]=r
    num = num +1
endwhile

close,lun
free_lun,lun
get_lun,lun
openr,lun,neutralfile
vec=dblarr(7)
while not eof(lun) do begin
    readf,lun,vec
    nh=vec[0]
    i=vec[1]
    j=vec[2]
    k=vec[3]
    theta=vec[4]
    phi=vec[5]
    r=vec[6]
    snhopher.density[i,j,k]=nh
endwhile
;    if sdataopher.phi[i,j,k] ne phi then stop & if sdataopher.theta[i,j,k] ne theta then stop & if sdataopher.r[i,j,k] ne r then stop ; in previous code replaced theta,phi, r
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
if modeler eq 'Heerikhuisen' then begin
    common grid_heerikhuisen,gridheer
    common nhdata_heerikhuisen,nhheer,snhheer 
    common data_heerikhuisen,dataheer,sdataheer
    nr=gridheer.r.num   ; number of r bins
    data = sdataheer
    data_neutral = snhheer
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

nn=0D
np=0D
tp=0D
theta=0D
phi=0D
r=0D
ur=0D
ut=0D ; Ut and Un used interchangably   

r_xyz=dblarr(ntheta+1,nphi,nr)
theta_xyz=dblarr(ntheta+1,nphi,nr)
phi_xyz=dblarr(ntheta+1,nphi,nr)

r_xyz[0:ntheta-1,0:nphi-1,0:nr-1]=data.r
theta_xyz[0:ntheta-1,0:nphi-1,0:nr-1]=data.theta
phi_xyz[0:ntheta-1,0:nphi-1,0:nr-1]=data.phi

r_xyz[ntheta,0:nphi-1,0:nr-1]=data.r[ntheta-1,0:nphi-1,0:nr-1]
theta_xyz[ntheta,0:nphi-1,0:nr-1]=180.
phi_xyz[ntheta,0:nphi-1,0:nr-1]=data.phi[ntheta-1,0:nphi-1,0:nr-1]

;x0=data.r*cos(data.phi*!pi/180.)*sin(data.theta*!pi/180.) 
;y0=data.r*sin(data.phi*!pi/180.)*sin(data.theta*!pi/180.)
;z0=data.r*cos(data.theta*!pi/180.)

x0=r_xyz*cos(phi_xyz*!pi/180.)*sin(theta_xyz*!pi/180.)
y0=r_xyz*sin(phi_xyz*!pi/180.)*sin(theta_xyz*!pi/180.)
z0=r_xyz*cos(theta_xyz*!pi/180.)

; Convert to ECL coords
x1=0.263023*x0-0.957015*y0+0.122333*z0
y1=0.96478*x0+0.261166*y0+0.0310853*z0
z1=-0.00227776*x0+0.126148*y0+0.992011*z0
r1=sqrt(x1^2.+y1^2.+z1^2.)
theta1=atan(y1/(x1+1e-10))
phi1=acos(z1/(r1+1e-10))

;-------------------
;this section for radial_ena_prod.pro
common rad_profile, radial_profile
rad_arr_3 = make_array(3,ntheta,nphi,nr,/double)
rad_arr = make_array(ntheta,nphi,nr,/double)
;;;;;;;;;;;;;;;;;;;;;;
fds = make_array(ntheta,nphi,/double)
;;;;;;;;;;;;;;;;;;;;;;
radial_profile = {radial_stuff, ena_flux:rad_arr_3, Temp:rad_arr_3, density:rad_arr_3, Ur:rad_arr, Ut:rad_arr, nH:rad_arr, emag:rad_arr, eram:rad_arr, vtheta:rad_arr, vphi:rad_arr}
gam = 5./3. ; adiabtic index
;--------------------
rlower=(r_lower-r_inner)/dr
rupper=(r_upper-r_inner)/dr

;lc=make_array(ntheta,nphi)
;cool=make_array(ntheta,nphi,nr,/double)
;vx=make_array(ntheta,nphi,nr,/double)
;vy=make_array(ntheta,nphi,nr,/double)
;vz=make_array(ntheta,nphi,nr,/double)
;read_stream, energy, lc, cool, vx, vy, vz

vena=energy_velocity(energy)
sig_cx=xsection(energy)
ext=cool*vena*sig_cx
ext=exp(ext)

;cool_length, cool

if extinction ne 1 then ext(*,*,*)=1.

; Creating wedge values
; wedge is region where plasmoid affects lobes
; Finding edges outside of wedge so can use edges to interpolate over wedge for given r
; wedge=1 means north of equator
; wedge=2 means south of equator
; Need to divide wedge into north/south based on collimation

boundaries_arr1=make_array(ntheta,nphi,nr,/double)
vas=make_array(ntheta,nphi,nr,/double)
psi_arr=make_array(ntheta,nphi,/double)
fi=make_array(ntheta,nphi,/double)
ni=make_array(ntheta,nphi,/double)
ni2=make_array(ntheta,nphi,/double)
ni_rat=make_array(ntheta,nphi,/double)
ni_rat2=make_array(ntheta,nphi,nr,/double)
ni_rat_nor=make_array(ntheta,nphi,/double)
ni_rat_plas=make_array(ntheta,nphi,/double)
np_ts=make_array(ntheta,nphi,/double)
r_ts=make_array(ntheta,nphi,/double)
np_psi=make_array(ntheta,nphi,/double)
r_arr=make_array(ntheta,nphi,/double)
n1_ts=make_array(ntheta,nphi,/double)
n2_ts=make_array(ntheta,nphi,/double)
nsw_ts=make_array(ntheta,nphi,/double)

bydi=make_array(ntheta,nphi,nr,/double)
;read_stream_japan, bydi

if vasyliunas eq 1 then begin ; to use vasyliunas PUI distribution
	for t=0,ntheta-1 do begin
		for p=0,nphi-1 do begin
			for rad=rlower,rupper do begin
				r=round(data.r[t,p,rad])
				theta=data.theta[t,p,rad]
				phi=data.phi[t,p,rad]
				ur=data.vr[t,p,rad]*1.e5
                                ut=data.vt[t,p,rad]*1.e5
				n_ism=0.1 ; ISM density at infinity (0.18)
;				n_ism=0.18
;				n_ism=0.015 ; Helium
;				n_ism=data_neutral.density[t,p,rad]
                                v0=20.e5 ; speed of Sun relative to neutral gas (25 km/s)
;				v0=25.e5
                                G=6.6726e-8
                                M_sun=1.989e33
                                tau_1au=1e-6 ; ionization rate at 1 AU = 8e-8 s^-1
;				tau_1au=4.7e-8 ; Helium
                                a_ref=au_cm ; reference distance where ionization rate is calculated for
                                lam=a_ref^2.*tau_1au/v0 ; characteristic distance
				mu=0.0
				if r eq 30 then begin
					r=2.
					for rinner=0,13 do begin ; calculating flux density within inner boundary
;					for rinner=0,3 do begin ; stops after r=8 AU, for r/lam~1 (for v0=20 km/s) 
						t_rad=theta*!pi/180.
						p_rad=phi*!pi/180.
						t_nose=data.theta[ntheta/2.,nphi/2.,0]*!pi/180. ; reference theta at nose
						p_nose=data.phi[ntheta/2.,nphi/2.,0]*!pi/180. ; reference phi at nose
						r_au=r*au_cm
						dr_au=dr*au_cm
						psi=acos(cos(t_rad)*cos(t_nose)+sin(t_rad)*sin(t_nose)*cos(abs(p_rad-p_nose))) ; angle from nose
	
						b1=sqrt((0.5*r_au*sin(psi))^2.+(1-mu)*G*M_sun/(v0^2.)*r_au*(1.-cos(psi)))+0.5*r_au*sin(psi)
						b2=sqrt((0.5*r_au*sin(psi))^2.+(1-mu)*G*M_sun/(v0^2.)*r_au*(1.-cos(psi)))-0.5*r_au*sin(psi)
						
						; analytically calculating db1 and db2
	                                        con=G*M_sun/(v0^2.)
        	                                db=con*(1-cos(psi))+0.5*r_au*(sin(psi))^2.
                	                        db=db/(2*sqrt(con*r_au*(1-cos(psi))+0.25*(r_au*sin(psi))^2.))
                        	                db1=db+0.5*sin(psi)
                                	        db2=db-0.5*sin(psi)
	
        	                                vas_1=db1*exp(-lam*psi/b1)
                	                        vas_2=db2*exp(-lam*(2.*!pi-psi)/b2)
                        	                vas[t,p,rad]=n_ism/(sin(psi))*(vas_1+vas_2)
                                	        psi_arr[t,p]=psi
	
        	                                fi[t,p]=tau_1au*(a_ref/r_au)^2.*dr_au*vas[t,p,rad]+fi[t,p]
						ni2[t,p]=fi[t,p]/3.e7
						r=r+dr				
					endfor
				endif
				if r lt 150 then begin
					if boundaries_arr[t,p,rad] eq 4 then begin
						t_rad=theta*!pi/180.
        					p_rad=phi*!pi/180.
        					t_nose=data.theta[ntheta/2.,nphi/2.,0]*!pi/180. ; reference theta at nose
        					p_nose=data.phi[ntheta/2.,nphi/2.,0]*!pi/180. ; reference phi at nose
        					r_au=r*au_cm
						dr_au=dr*au_cm
        					psi=acos(cos(t_rad)*cos(t_nose)+sin(t_rad)*sin(t_nose)*cos(abs(p_rad-p_nose))) ; angle from nose 

	        				b1=sqrt((0.5*r_au*sin(psi))^2.+G*M_sun/(v0^2.)*r_au*(1.-cos(psi)))+0.5*r_au*sin(psi)
        					b2=sqrt((0.5*r_au*sin(psi))^2.+G*M_sun/(v0^2.)*r_au*(1.-cos(psi)))-0.5*r_au*sin(psi)

						; analytically calculating db1 and db2
						con=G*M_sun/(v0^2.)
						db=con*(1-cos(psi))+0.5*r_au*(sin(psi))^2.
						db=db/(2*sqrt(con*r_au*(1-cos(psi))+0.25*(r_au*sin(psi))^2.))
						db1=db+0.5*sin(psi)
						db2=db-0.5*sin(psi)

        					vas_1=db1*exp(-lam*psi/b1)
        					vas_2=db2*exp(-lam*(2.*!pi-psi)/b2)
 	       					vas[t,p,rad]=n_ism/(sin(psi))*(vas_1+vas_2)
        					psi_arr[t,p]=psi
						r_arr[t,p]=r
						nh_avg_ct_nose=where(boundaries_arr[ntheta/2.,nphi/2.,0:60] eq 4)
						nh_avg_nose=mean(data_neutral.density[ntheta/2.,nphi/2.,nh_avg_ct_nose])*2.
        					fi[t,p]=tau_1au*(a_ref/r_au)^2.*dr_au*vas[t,p,rad]+fi[t,p]
; 	      					fi[t,p]=tau_1au*(a_ref/r_au)^2.*dr_au*nn+fi[t,p]
;       	 				ni[t,p]=fi[t,p]/sqrt(ur^2.+ut^2.)
						ni[t,p]=fi[t,p]/3.e7 ; 3.e7 = mean of speed just before TS
						if r eq 78 then np_psi[t,p]=data.density[t,p,rad]
						if boundaries_arr[t,p,rad] ne 3 and boundaries_arr[t,p,rad+1] eq 3 then r_ts[t,p]=data.r[t,p,rad+1]
						if boundaries_arr[t,p,rad] ne 3 and boundaries_arr[t,p,rad+1] ne 3 and boundaries_arr[t,p,rad+2] eq 3 then r_ts[t,p]=data.r[t,p,rad+2]
						if boundaries_arr[t,p,rad] ne 3 and boundaries_arr[t,p,rad+2] ne 3 and boundaries_arr[t,p,rad+3] eq 3 then r_ts[t,p]=data.r[t,p,rad+3] 
	       				endif
				endif
			endfor
		endfor
	endfor 				
	ni[ntheta/2.,nphi/2.]=(ni[ntheta/2.-1.,nphi/2.]+ni[ntheta/2.+1.,nphi/2.]+ni[ntheta/2.,nphi/2.-1.]+ni[ntheta/2.,nphi/2.+1.])/4. ; for where psi-->0
	ni[ntheta/2.,0]=(ni[ntheta/2.,1]+ni[ntheta/2.,nphi-1.]+ni[ntheta/2.-1.,0]+ni[ntheta/2.+1,0])/4. ; for where psi-->pi
	tail=where(psi_arr gt 2.90,tct)
	fixtail=where(psi_arr le 2.90 and psi_arr gt 2.84,fct)
	ni[tail]=ni[fixtail[0]]

        ni2[ntheta/2.,nphi/2.]=(ni2[ntheta/2.-1.,nphi/2.]+ni2[ntheta/2.+1.,nphi/2.]+ni2[ntheta/2.,nphi/2.-1.]+ni2[ntheta/2.,nphi/2.+1.])/4. ; for where psi-->0
        ni2[ntheta/2.,0]=(ni2[ntheta/2.,1]+ni2[ntheta/2.,nphi-1.]+ni2[ntheta/2.-1.,0]+ni2[ntheta/2.+1,0])/4. ; for where psi-->pi
        ni2[tail]=ni2[fixtail[0]]

endif

psi_arr=psi_arr*180/!pi
psi_arr=round(psi_arr)
n_av=fltarr(181)
r_av=fltarr(181)
for n=0,180 do begin
	w=where(psi_arr eq n)
	n_av[n]=mean(np_psi[w])
	r_av[n]=mean(r_arr[w])
endfor

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
;vas=make_array(ntheta,nphi,nr,/double)
;psi_arr=make_array(ntheta,nphi,/double)
;fi=make_array(ntheta,nphi,/double)
;ni=make_array(ntheta,nphi,/double)
;ni_rat=make_array(ntheta,nphi,/double)
psi_arr1=make_array(ntheta,nphi,nr,/double)
;if pui_model eq 'Malama' then a_z=make_array(ntheta,nphi,nr,/double)
stream_nose={nose,theta:fltarr(nr),phi:fltarr(nr),r:fltarr(nr),dens:fltarr(nr),temp:fltarr(nr),vr:fltarr(nr),vt:fltarr(nr),ext:fltarr(nr),nh:fltarr(nr),rat:fltarr(nr)}
stream_tail={tail,theta:fltarr(nr),phi:fltarr(nr),r:fltarr(nr),dens:fltarr(nr),temp:fltarr(nr),vr:fltarr(nr),vt:fltarr(nr),ext:fltarr(nr),nh:fltarr(nr),rat:fltarr(nr)}
stream_lobe={lobe,theta:fltarr(nr),phi:fltarr(nr),r:fltarr(nr),dens:fltarr(nr),temp:fltarr(nr),vr:fltarr(nr),vt:fltarr(nr),ext:fltarr(nr),nh:fltarr(nr),rat:fltarr(nr)}

if flag_survival eq 1 then ionization_rate,energy,bet ; reads in ionization rate
if flag_survival eq 0 then bet=dblarr(ntheta,nphi,nr)

;dens_file='malama_density_interpolated.txt' ; reading in PUI densities from Malama
;readcol, dens_file, z_np, np1_f, np2_f, np3_f, np4_f, /silent

;temp_file='malama_temperature_interpolated.txt' ; reading in PUI temps from Malama
;readcol, temp_file, z_tp, tp1_f, tp2_f, tp3_f, tp4_f, /silent

; Since Malama has TS at 112 AU, bringing TS to our model's TS
;np1_uw=dblarr(nr)
;np2_uw=dblarr(nr)
;np3_uw=dblarr(nr)
;np4_uw=dblarr(nr)
;tp1_uw=dblarr(nr)
;tp2_uw=dblarr(nr)
;tp3_uw=dblarr(nr)
;tp4_uw=dblarr(nr)
;for rad=rlower, rupper do begin 
;	np1_uw[rad]=np1_f[rad+10]
;	np2_uw[rad]=np2_f[rad+10]
;	np3_uw[rad]=np3_f[rad+10]
;	np4_uw[rad]=np4_f[rad+10]
;	tp1_uw[rad]=tp1_f[rad+10]
;	tp2_uw[rad]=tp2_f[rad+10]
;	tp3_uw[rad]=tp3_f[rad+10]
;	tp4_uw[rad]=tp4_f[rad+10]
;endfor

percentage=-1
for t=0, ntheta-1 do begin
    for p=0, nphi-1 do begin
	
	depth=where(abs(1/exp(1)-ext[t,p,*]) eq min(abs(1/exp(1)-ext[t,p,*])),dct) ; cooling depth calc
        lc[t,p]=min(depth)

	for rad=rlower,rupper do begin
;	for rad=rlower, lc[t,p] do begin
		if radial_range eq 0. then begin ; FOR TS TO HP - MZK
;                        if boundaries_arr[t,p,rad] eq 3 then begin; REGION 3 = HELIOSHEATH - MZK
			if boundaries_arr[t,p,rad] eq 3 or boundaries_arr[t,p,rad] eq 6 then begin; REGION 3 = HELIOSHEATH - MZK

;print, "Begin make_flux"
;print, "radial_range=", radial_range
;print, "theta=", t
;print, "phi=", p
;print, "rad=", rad
;print, "rupper=", rupper
;print, "rlower=", rlower
;print, "region:", boundaries_arr[t,p,rad] 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ASIDE FROM 2 TABS PER LINE, REST IS UNCHANGED UNTIL ENDIF STATEMENTS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		                if modeler eq 'Opher' then nn=.1
               		        if modeler eq 'Heerikhuisen' then nn=data_neutral.density[t,p,rad]
				if modeler eq 'Opher_withneutrals' or modeler eq 'Opher_multiion_withneutrals' then nn=data_neutral.density[t,p,rad]
             
		;-------------------
                ;;this section for radial_ena_prod.pro
                ;;FLOW SPEEDS CHANGED HERE
    		;if (data.r[t,p,rad] ge 130.) then begin  ; data.r is in AU 
                ;    ur = 0.
                ;    ut=0.;data.vt[t,p,rad]*10.^5.
		;    uphi = 0. ; Un in Voyager
                ;    utheta = 0.;-26.*10.^5. ;UT in Voyager
                ;endif else begin
                ;-------------------     
                  ;REGULAR FLOWS USED, do not comment out inside of this portion of loop
            		       ur = data.vr[t,p,rad]*10.^5. 
                  	       ut = data.vt[t,p,rad]*10.^5.
	       	               uphi = data.vphi[t,p,rad]*10.^5.  ; redunant to Ut, for when Merav wanted to split apart non radial component
                  	       utheta = data.vtheta[t,p,rad]*10.^5.; redunant to Ut, for when Merav wanted to split apart non radial component
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
flux[t,p]=calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0,bet[t,p,rad])+flux[t,p] ; kappa0 and bet by MZK  
radial_profile.ena_flux(0,t,p,rad) = calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0,bet[t,p,rad]) ; trying to find single ion flux, kappa0 and bet by MZK
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      	 
;				flux[t,p]=calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0)+flux[t,p] ; kappa0 by MZK  	
;				radial_profile.ena_flux(0,t,p,rad) = calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0) ; trying to find single ion flux, kappa0 by MZK

;print, "post-flux calculation"		
;print, "ur=",ur                                         
;print, "ut=",ut                                         
;print, "uphi=",uphi
;print, "utheta=",utheta                                 
;print, "theta (deg)=",theta                                     
;print, "phi(deg)=",phi                                  
;print, "radius (AU)=",r                                         
;print, "np=", np
;print, "tp=", tp         				
;print, "point flux=", radial_profile.ena_flux(0,t,p,rad)
;print, "total flux=", flux[t,p]
;print, "" 
                 		endif else begin ;MULTIION FLUID CASE ------------------------------------------------------------------------

                        		np=data.density[t,p,rad]
                			tp=data.temp[t,p,rad]
   					nsw=data.density_sw[t,p,rad]
; 					npu3=data.density_pu3[t,p,rad]
 					npt=data.density_pt[t,p,rad]
					npr=data.density_pr[t,p,rad]
                       ;-------------------
                       ;this section for radial_ena_prod.pro
                       ; be careful of parantheses when changing different runs
                       ;if (data.r[t,p,rad] gt 130.) and (data.r[t,p,rad] le r_hp) then begin   
                        ;  tsw = data.temp_sw[t,p,rad] + A_temp*(gamma-1.)*[.5*mp/kb*data.density[t,p,rad]/data.density_sw[t,p,rad]*(data.vr[t,p,rad]^2.*(10^5.)^2. + data.vt[t,p,rad]^2.*(10^5.)^2.)];+1.*(gamma-1.)*((ophermagnetic.Bx[t,p,rad]^2. +ophermagnetic.By[t,p,rad]^2. +ophermagnetic.Bz[t,p,rad]^2.)/(8.*!pi)/kb/data.density_sw[t,p,rad]) ;+A_temp*(gamma-1.)*([.5*mp/kb*data.density[t,p,rad]/data.density_sw[t,p,rad]*(data.vr[t,p,rad]^2.)*(10^5.)^2.]) ;+ data.vt[t,p,rad]^2.)*(10^5.)^2.]) ;+ B_temp*(gamma-1.)*((ophermagnetic.Bx[t,p,rad]^2. +ophermagnetic.By[t,p,rad]^2. +ophermagnetic.Bz[t,p,rad]^2.)/(8.*!pi)/kb/data.density_sw[t,p,rad]) 
                         ; tpu3 = data.temp_pu3[t,p,rad] + C_temp*(gamma-1.)*[.5*mp/kb*data.density[t,p,rad]/data.density_pu3[t,p,rad]*(data.vr[t,p,rad]^2.*(10^5.)^2. + data.vt[t,p,rad]^2.*(10^5.)^2.)] ;+ D_temp*(gamma-1.)*((ophermagnetic.Bx[t,p,rad]^2. +ophermagnetic.By[t,p,rad]^2. +ophermagnetic.Bz[t,p,rad]^2.)/(8.*!pi)/kb/data.density_pu3[t,p,rad])
                       ;-------------------
                       ;endif else begin
                          ;TEMPERATURES FROM MODEL RESULTS USED, do not comment out inside of this portion of loop
                          		tsw  =  data.temp_sw[t,p,rad] 
;                          		tpu3 =  data.temp_pu3[t,p,rad]
					tpt = data.temp_pt[t,p,rad]
					tpr = data.temp_pr[t,p,rad]
					if boundaries_arr[t,p,rad] eq 3. then kappa0=kappa 
					if boundaries_arr[t,p,rad] eq 6. then kappa0=kappaOH
                      ;endelse 

                        ;This is used for making countrate files. For regular maps make sure flows and temperatures are just as pulled from data file. Do not change.  
;					flux[t,p,0]=calc_diff_flux_multiion(dist,energy,r,tp,ur,ut,nn,np,theta,Uphi,Utheta,kappa0)+flux[t,p,0] ; total ion, kappa0 by MZK
;                        		flux[t, p,1]=calc_diff_flux_multiion(dist,energy,r,tsw,ur,ut,nn,nsw,theta,Uphi,Utheta,kappa0)+flux[t,p,1] ; SW component, kappa0 by MZK
;                        		flux[t,p,2]= calc_diff_flux_multiion(dist,energy,r,tpu3,ur,ut,nn,npu3,theta,Uphi,Utheta,kappa0)+flux[t,p,2] ; PUI component, kappa0 by MZK
					flux[t,p,0]=calc_diff_flux(dist,energy,r,tsw,ur,ut,nn,nsw,theta,kappa0,bet[t,p,rad])+flux[t,p,0] ; total ion, kappa0 by MZK
                                        flux[t,p,1]=calc_diff_flux(dist,energy,r,tpt,ur,ut,nn,npt,theta,kappa0,bet[t,p,rad])+flux[t,p,1] ; SW component, kappa0 by MZK
                                        flux[t,p,2]= calc_diff_flux(dist,energy,r,tpr,ur,ut,nn,npr,theta,kappa0,bet[t,p,rad])+flux[t,p,2] ; PUI component, kappa0 by MZK
;                      			flux[t,p,3]=flux[t,p,0]+flux[t,p,1]+flux[t,p,2]  
                       ;-------------------
                       ;this section for radial_ena_prod.pro
                       ;For changing the distribution function for select regions, takes advantage of passed variables A_temp and B_temp
                       ;if (data.r[t,p,rad] ge r_hp) and (data.r[t,p,rad] le 350.) then begin 
                       ; kappa = B_temp 
		       ; dist = 'Kappa'+string(kappa,format='(f4.2)')
	               ;endif else dist  = 'Maxwellian'
                        ;ALWAYS LEAVE radial_profile.ena_flux uncommented - there are 3

                       	 		radial_profile.ena_flux(0,t,p,rad) = calc_diff_flux(dist,energy,r,tsw,ur,ut,nn,nsw,theta,kappa0,bet[t,p,rad]) ; kappa0 by MZK
                        		radial_profile.ena_flux(1,t,p,rad) = calc_diff_flux(dist,energy,r,tpt,ur,ut,nn,npt,theta,kappa0,bet[t,p,rad]) ; kappa0 by MZK

                       ; if (data.r[t,p,rad] le r_hp) then begin 
                       ; kappa = A_temp 
		       ; dist = 'Kappa'+string(kappa,format='(f4.2)')
	               ;endif else dist  = 'Maxwellian'
 		        		radial_profile.ena_flux(2,t,p,rad) = calc_diff_flux(dist,energy,r,tpr,ur,ut,nn,npr,theta,kappa0,bet[t,p,rad]) ; kappa0 by MZK
                       ;When done best to set to default, as it is a passed global parameter
                       ;dist = 'Maxwellian'
                       
                        ;Track parameters for globally passed variable
  ;                      		radial_profile.Temp(0,t,p,rad) = tp & radial_profile.Temp(1,t,p,rad) = tsw &radial_profile.Temp(2,t,p,rad) = tpu3
 ;           				radial_profile.density(0,t,p,rad) = np & radial_profile.density(1,t,p,rad) = nsw &radial_profile.density(2,t,p,rad) = npu3
					radial_profile.nH(t,p,rad) = nn  & radial_profile.Ur(t,p,rad) =data.vr[t,p,rad]^2.  & radial_profile.Ut(t,p,rad) =  data.vt[t,p,rad]
                        		radial_profile.vtheta(t,p,rad) = utheta  & radial_profile.vphi(t,p,rad) = uphi
					radial_profile.eram(t,p,rad ) = [.5*mp*data.density[t,p,rad]*(data.vr[t,p,rad]^2. + data.vt[t,p,rad]^2.)*(10^5.)^2.]
           				radial_profile.emag(t,p,rad) = (ophermagnetic.Bx[t,p,rad]^2. +ophermagnetic.By[t,p,rad]^2. +ophermagnetic.Bz[t,p,rad]^2.)/(8.*!pi)
                        ;------------------- 
				endelse
			endif ; ADDED BY MZK FOR REGION IF STATEMENT
		endif ; ADDED BY MZK FOR RADIAL RANGE IF STATEMENT

;\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
;/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    RADIAL RANGE 1      \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
;\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

; THE FOLLOWING HAS BEEN ADDED BY MZK TO ALLOW FOR THE RADIAL RANGE TO DICTATE THE ENA MAP RANGES
;

                if radial_range eq 1. then begin ; FOR HP to 1500 AU - MZK
                        if boundaries_arr[t,p,rad] eq 2 or boundaries_arr[t,p,rad] eq 1 then begin ; REGION 2 = BEYOND HP - MZK

                                if modeler eq 'Opher' then nn=.1
                                if modeler eq 'Heerikhuisen' then nn=data_neutral.density[t,p,rad]
                                if modeler eq 'Opher_withneutrals' or modeler eq 'Opher_multiion_withneutrals' then nn=data_neutral.density[t,p,rad]


                  ;REGULAR FLOWS USED, do not comment out inside of this portion of loop
                               ur = data.vr[t,p,rad]*10.^5.
                               ut = data.vt[t,p,rad]*10.^5.
                               uphi = data.vphi[t,p,rad]*10.^5.  ; redunant to Ut, for when Merav wanted to split apart non radial component
                               utheta = data.vtheta[t,p,rad]*10.^5.; redunant to Ut, for when Merav wanted to split apart non radial component
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
flux[t,p]=calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0,bet[t,p,rad])+flux[t,p] ; kappa0 and bet by MZK  
radial_profile.ena_flux(0,t,p,rad) = calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0,bet[t,p,rad]) ; trying to find single ion flux, kappa0 and bet by MZK
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;            
;                                        flux[t,p]=calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0)+flux[t,p]; kappa0 by MZK
;					radial_profile.ena_flux(0,t,p,rad) = calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0) ; trying to find single ion flux, kappa0 by MZK			

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

                        ;This is used for making countrate files. For regular maps make sure flows and temperatures are just as pulled from data file. Do not change.  
                                        flux[t,p,0]=calc_diff_flux_multiion(dist,energy,r,tp,ur,ut,nn,np,theta,Uphi,Utheta,kappa0)+flux[t,p,0] ; total ion, kappa0 by MZK
                                        flux[t,p,1]=calc_diff_flux_multiion(dist,energy,r,tsw,ur,ut,nn,nsw,theta,Uphi,Utheta,kappa0)+flux[t,p,1] ; SW component, kappa0 by MZK
                                        flux[t,p,2]= calc_diff_flux_multiion(dist,energy,r,tpu3,ur,ut,nn,npu3,theta,Uphi,Utheta,kappa0)+flux[t,p,2] ; PUI component, kappa0 by MZK

                        ;ALWAYS LEAVE radial_profile.ena_flux uncommented - there are 3
                                        radial_profile.ena_flux(0,t,p,rad) = calc_diff_flux_multiion(dist,energy,r,tp,ur,ut,nn,np,theta,Uphi,Utheta,kappa0) ; kappa0 by MZK
                                        radial_profile.ena_flux(1,t,p,rad) = calc_diff_flux_multiion(dist,energy,r,tsw,ur,ut,nn,nsw,theta,Uphi,Utheta,kappa0) ; kappa0 by MZK
                                        radial_profile.ena_flux(2,t,p,rad) = calc_diff_flux_multiion(dist,energy,r,tpu3,ur,ut,nn,npu3,theta,Uphi,Utheta,kappa0) ; kappa0 by MZK

                        ;Track parameters for globally passed variable
                                        radial_profile.Temp(0,t,p,rad) = tp & radial_profile.Temp(1,t,p,rad) = tsw &radial_profile.Temp(2,t,p,rad) = tpu3
                                        radial_profile.density(0,t,p,rad) = np & radial_profile.density(1,t,p,rad) = nsw &radial_profile.density(2,t,p,rad) = npu3
                                        radial_profile.nH(t,p,rad) = nn  & radial_profile.Ur(t,p,rad) =data.vr[t,p,rad]^2.  & radial_profile.Ut(t,p,rad) =  data.vt[t,p,rad]
                                        radial_profile.vtheta(t,p,rad) = utheta  & radial_profile.vphi(t,p,rad) = uphi
                                        radial_profile.eram(t,p,rad ) = [.5*mp*data.density[t,p,rad]*(data.vr[t,p,rad]^2. + data.vt[t,p,rad]^2.)*(10^5.)^2.]
                                        radial_profile.emag(t,p,rad) = (ophermagnetic.Bx[t,p,rad]^2. +ophermagnetic.By[t,p,rad]^2. +ophermagnetic.Bz[t,p,rad]^2.)/(8.*!pi)
                        ;------------------- 
                                endelse
			endif
		endif

;<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
;*************************************************       RADIAL RANGE 2       *************************************************
;<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

;		phi=data.phi[t,p,rad]
                if radial_range eq 2. then begin ; FOR TS to 1500 AU - MZK
			if boundaries_arr[t,p,rad] eq 3 or boundaries_arr[t,p,rad] eq 2 or boundaries_arr[t,p,rad] eq 1 or boundaries_arr[t,p,rad] eq 6 or boundaries_arr[t,p,rad] eq 4 then begin ; REGION 3 = TS TO BEYOND- MZK

			; added to restrict maps to downstream only
;			if (phi le 90. or phi ge 270.) and (boundaries_arr[t,p,rad] eq 3 or boundaries_arr[t,p,rad] eq 2 or boundaries_arr[t,p,rad] eq 1) then begin

;tp =  data.temp[t,p,rad] ; added to restrict region 3 to ISM between lobes
;                        if boundaries_arr[t,p,rad] eq 3 and tp gt 2.e5 or boundaries_arr[t,p,rad] eq 2 and tp gt 2.e5 or boundaries_arr[t,p,rad] eq 1 and tp gt 2.e5 then begin ; FOR HS and ISM between lobes

                                if modeler eq 'Opher' then nn=.1
                                if modeler eq 'Heerikhuisen' then nn=data_neutral.density[t,p,rad]
                                if modeler eq 'Opher_withneutrals' or modeler eq 'Opher_multiion_withneutrals' then nn=data_neutral.density[t,p,rad]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; For modification with neutrals as separate populations - MZK
;
;				if modeler eq 'Opher_withneutrals' then begin
;					nn1=data_neutral.density[t,p,rad,0]
;					nn2=data_neutral.density[t,p,rad,1]
;					nn3=data_neutral.density[t,p,rad,2]
;					nn4=data_neutral.density[t,p,rad,3]
;				endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  ;REGULAR FLOWS USED, do not comment out inside of this portion of loop
                               ur = data.vr[t,p,rad]*10.^5.
                               ut = data.vt[t,p,rad]*10.^5.
                               uphi = data.vphi[t,p,rad]*10.^5.  ; redunant to Ut, for when Merav wanted to split apart non radial component
                               utheta = data.vtheta[t,p,rad]*10.^5.; redunant to Ut, for when Merav wanted to split apart non radial component
                ;endelse 

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
radial_profile.ena_flux(0,t,p,rad) = calc_diff_flux(dist,energy,r,tp,ur,ut,nn,np,theta,kappa0,bet[t,p,rad]); trying to find single ion flux, kappa0 and bet by MZK


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
;;;;;;;;
;
; Using Malama PUI populations
; PUI 1: Created in supersonic SW from primary/secondary interstellar H
; PUI 2: Created in supersonic SW from ENAs
; PUI 3: Created in IHS from primary/secondary interstellar H
; PUI 4: Created in IHS from ENAs

if boundaries_arr[t,p,rad] eq 3 and r le 96. then begin
	np_rat[t,p]=np/data.density[15,33,rad]
	tp_rat[t,p]=tp/data.temp[15,33,rad]
endif

np_scl=data.density[15,33,31]/0.00324
tp_scl=data.temp[15,33,31]/3.133e6

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

if vasyliunas eq 1 then ni_rat[t,p]=(nh_avg[t,p]/nh_avg_nose)*(r_ts[t,p]/r_ts[ntheta/2.,nphi/2.]) ; using nH and r -- Kornbleuth et al. (2018)
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

        nsw_frac=1-ntr_frac-nref_frac;-nen_frac
        ;esw_frac=1-etr_frac-eref_frac;-een_frac ; turn off when using Zank method

        ttr_frac=etr_frac/ntr_frac
        tref_frac=eref_frac/(nref_frac*(1-nen_frac))
        ;ten_frac=1e-15;een_frac/nen_frac ; turn off when using Zank method
        tsw_frac=0.04/nsw_frac;esw_frac/nsw_frac
        ten_frac=(1-etr_frac-eref_frac-0.04)/(nen_frac*nref_frac)
        if colde eq 1 and giacalone eq 1 then tp=2*tp/(1+tsw_frac)
endif

bsw=sqrt(ophermagnetic.Bx[t,p,rad]^2.+ophermagnetic.By[t,p,rad]^2.+ophermagnetic.Bz[t,p,rad]^2.)
;if bsw gt 0.20 then np=5.*np

if t eq 5 or t eq 6 then if p ge 42 and p le 46 then if boundaries_arr[t,p,rad] eq 3 and boundaries_arr[t,p,rad+1] eq 2 then np=np*20.
if t eq 5 or t eq 6 then if p ge 42 and p le 46 then if boundaries_arr[t,p,rad] eq 3 and boundaries_arr[t,p,rad+1] eq 2 then tp=tp*10.

if finite(ni_rat2[t,p,rad]) ne 1 then stop
zirn_rat[t,p]=a[t,p]/a[ntheta/2.,nphi/2.]
	if pui_model eq 'Malama' then begin
                ; Primary Method for Modeling PUIs
                if vasyliunas eq 1 then np1=ni_rat2[t,p,rad]*ntr_frac*np*ext[t,p,rad] else np1=ntr_frac*np*ext[t,p,rad] ; MHD Profile w/ Vasyliunas scaling
		if vasyliunas eq 1 then np2=ni_rat2[t,p,rad]*nref_frac*np*ext[t,p,rad] else np2=nref_frac*(1-nen_frac)*np*ext[t,p,rad] ; MHD Profile w/ Vasyliunas scaling
		if giacalone ne 1 then nsw=(1-ni_rat2[t,p,rad]*(ntr_frac+nref_frac))*np*ext[t,p,rad] else nsw=nsw_frac*np*ext[t,p,rad] ; solar wind
                if giacalone ne 1 then np3=np-nsw-np1-np2 else np3=nref_frac*nen_frac*np*ext[t,p,rad]; np3=nen_frac*np*ext[t,p,rad]; Locally created PUIs
                if np3 le 0 then np3=1e-10
		;if boundaries_arr[t,p,rad] eq 6 and t gt 1 then np3=np3 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then np3=1e-10
		;if boundaries_arr[t,p,rad] eq 6 and t gt 1 then np1=np1 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then np1=1e-10
		;if boundaries_arr[t,p,rad] eq 6 and t gt 1 then np2=np2 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then np2=1e-10
		;if boundaries_arr[t,p,rad] eq 6 and t gt 1 then nsw=nsw else if boundaries_arr[t,p,rad] eq 6 and t le 1 then nsw=1e-10
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
		;if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tp3=(bsw/0.1)^2.*tp3 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tp3=1e-10;3.*tp
		;if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tp1=(bsw/0.1)^2.*tp1 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tp1=1e-10; 1e-10
		;if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tp2=(bsw/0.1)^2.*tp2 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tp2=1e-10; 1e-10
		;if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tsw=(bsw/0.1)^2.*tsw else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tsw=1e-10; 1e-10
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tp3=tp3 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tp3=1e-10;3.*tp
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tp1=tp1 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tp1=1e-10; 1e-10
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tp2=tp2 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tp2=1e-10; 1e-10
                if boundaries_arr[t,p,rad] eq 6 and t gt 1 then tsw=(1+(bsw/0.15)^2.)*tp1 else if boundaries_arr[t,p,rad] eq 6 and t le 1 then tsw=1e-10; 1e-10

	endif else if pui_model eq 'Zirnstein' then begin
		if a_z[t,p,rad] eq 0 and boundaries_arr[t,p,rad] eq 3 then a_z[t,p,rad]=1e-10
 
		; Zirnstein densities
		np=data.density[t,p,rad]

                nsw=(1-a_z[t,p,rad])*np*ext[t,p,rad] ; Solar wind
                np1=a_z[t,p,rad]*(1-b)*np*ext[t,p,rad] ; Transmitted PUIs
                np2=a_z[t,p,rad]*b*np*ext[t,p,rad] ; Reflected PUIs

		if np-nsw-np1-np2 ge 0 then np3=np-nsw-np1-np2 else np3=1e-10 ; Locally created PUIs
		if nsw eq 0 then nsw=1e-10
;		if np1 eq 0 then np1=1e-10
;		if np2 eq 0 then np2=1e-10

		; Zirnstein temperatures
		tp=data.temp[t,p,rad];*0.8
;		tp=2*tp/(1+0.04/(1-a[t,p])) ; alternate method for defining temps from Zirnstein et al. (2017)
;		tsw=0.04*tp/(1-a[t,p]) ; Solar wind
;		tp1=0.5*tp/(a[t,p]*(1-b)) ; Transmitted PUIs
;		tp2=0.46*tp/(a[t,p]*b) ; Reflected PUIs
;		tp3=tsw ; Locally created PUIs

                if colde eq 1 then tp=2.*tp/(1+esw_frac/(1-a_z[t,p,rad])) else tp=tp/2.; alternate method for defining temps from Zirnstein et al. (2017)
                tsw=esw_frac*tp/(1-a_z[t,p,rad]) ; Solar wind
                tp1=etr_frac*tp/(a_z[t,p,rad]*(1-b)) ; Transmitted PUIs
                tp2=eref_frac*tp/(a_z[t,p,rad]*b) ; Reflected PUIs
                tp3=tsw ; Locally created PUIs
;                tp=2*tp/(1+0.04/(1-a_z[t,p,rad])) ; alternate method for defining temps from Zirnstein et al. (2017)
;		tsw=0.04*tp/(1-a_z[t,p,rad])
;		tp1=0.36*tp/(a_z[t,p,rad]*(1-b))
;		tp2=0.60*tp/(a_z[t,p,rad]*b)
;		tp3=tsw

		
		if boundaries_arr[t,p,rad-1] ne 3 then tp_ts[t,p]=tp1
	endif else if pui_model eq 'Other' then begin
		; Test case
		np=data.density[7,5,rad]
		nsw=np
		np1=np
		np2=np
		np3=np

		; Test case
		tp=data.temp[7,5,rad]
		tsw=snhopher.temp[t,p,rad]
		tp1=tp*4.
		tp2=snhopher.temp[t,p,rad]
		tp3=snhopher.temp[t,p,rad]
	endif

if tp1 le 0 or np1 le 0 or tp2 le 0 or np2 le 0 or tp3 le 0 or np3 le 0 or tsw le 0 or nsw le 0 then stop
if tp1 ne tp1 or np1 ne np1 or tp2 ne tp2 or np2 ne np2 or tp3 ne tp3 or np3 ne np3 or tsw ne tsw or nsw ne nsw then stop

endif else begin 
        np1=1.e-15
        np2=1.e-15
        np3=1.e-15
        nsw=1.e-15

        tp1=1.
        tp2=1.
        tp3=1.
        tsw=1.

	if giacalone eq 1 then een_frac=0.
	if giacalone eq 1 then em_frac=0.
	if giacalone eq 1 then vshock=0.
endelse

; For testing constant value in heliosphere for np1, tp1, u
;np1=4e-4
;np1=np1*1.22
;nn=0.1
;tp1=1e7
;ur=100e5
;ut=0.

dens_test1[t,p,rad]=np1
temp_test1[t,p,rad]=tp1
dens_test2[t,p,rad]=np2
temp_test2[t,p,rad]=tp2
dens_test3[t,p,rad]=np3
temp_test3[t,p,rad]=tp3
dens_test_sw[t,p,rad]=nsw
temp_test_sw[t,p,rad]=tsw
temp_test_p[t,p,rad]=tp
psi_arr1[t,p,rad]=psi_arr[t,p]
if vas[t,p,rad] ne vas[t,p,rad] then vas[t,p,rad]=0.

data.vr[t,p,rad]=ur
data.vt[t,p,rad]=ut

v_plasma = plasma_frame(energy,ur,ut)
f_test1[t,p,rad]=double(maxwellian_fn(tp1,np1,v_plasma))
v_plas[t,p,rad]=v_plasma
if rad eq lc[t,p] then lc_test[t,p,rad]=1. else lc_test[t,p,rad]=0
if boundaries_arr[t,p,rad] eq 3 and boundaries_arr[t,p,rad-1] eq 4 then shock_strength[t,p]=up[t,p,rad-2]/up[t,p,rad]

if t eq 15 and p eq 30 and boundaries_arr[t,p,rad] eq 3 then begin
	stream_nose.theta[rad]=theta
	stream_nose.phi[rad]=phi
	stream_nose.r[rad]=r
	stream_nose.dens[rad]=np
	stream_nose.temp[rad]=data.temp[t,p,rad]/2.
	stream_nose.vr[rad]=ur
	stream_nose.vt[rad]=ut
	stream_nose.ext[rad]=ext[t,p,rad]
	stream_nose.nh[rad]=nn
	stream_nose.rat[rad]=ni_rat2[t,p,rad]
endif

if t eq 15 and p eq 00 and boundaries_arr[t,p,rad] eq 3 then begin
        stream_tail.theta[rad]=theta
        stream_tail.phi[rad]=phi
        stream_tail.r[rad]=r
        stream_tail.dens[rad]=np
        stream_tail.temp[rad]=data.temp[t,p,rad]/2.
        stream_tail.vr[rad]=ur
        stream_tail.vt[rad]=ut
        stream_tail.ext[rad]=ext[t,p,rad]
        stream_tail.nh[rad]=nn
        stream_tail.rat[rad]=ni_rat2[t,p,rad]
endif

if t eq 10 and p eq 0 and boundaries_arr[t,p,rad] eq 3 then begin
        stream_lobe.theta[rad]=theta
        stream_lobe.phi[rad]=phi
        stream_lobe.r[rad]=r
        stream_lobe.dens[rad]=np
        stream_lobe.temp[rad]=data.temp[t,p,rad]/2.
        stream_lobe.vr[rad]=ur
        stream_lobe.vt[rad]=ut
        stream_lobe.ext[rad]=ext[t,p,rad]
        stream_lobe.nh[rad]=nn
        stream_lobe.rat[rad]=ni_rat2[t,p,rad]
endif

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
; Tests for PUI density ratio, mass flux, etc.
;                                        if boundaries_arr[t,p,rad] eq 3 and boundaries_arr[t,p,rad-1] eq 4 then flux[t,p,0]=a_z[t,p,rad]
;					if boundaries_arr[t,p,rad] eq 3 and boundaries_arr[t,p,rad-26] eq 4 and boundaries_arr[t,p,rad-25] eq 3 then flux[t,p,0]=np*sqrt(ur^2+ut^2)/1e5
;					if flux[t,p,0] ne 0 and p gt 14 and p lt 45 then flux[t,p,0]=0
;					if boundaries_arr[t,p,rad] eq 3 and boundaries_arr[t,p,rad-1] eq 3 and boundaries_arr[t,p,rad-2] eq 4 then flux[t,p,0]=data.density[t,p,rad]/data.density[t,p,rad-4]
;                                       if boundaries_arr[t,p,rad] eq 3 and boundaries_arr[t,p,rad-1] eq 3 and boundaries_arr[t,p,rad-2] eq 4 then flux[t,p,0]=data.density[t,p,rad]/data.density[t,p,rad-4]

              				; Test for Constant ENA Production
;					if sqrt(x0[t,p,rad]^2.+y0[t,p,rad]^2.) le 111. and data.r[t,p,rad] gt 66 then boundaries_arr1[t,p,rad] = 3 $
;					else if data.r[t,p,rad] le 66. then boundaries_arr1[t,p,rad]=4 else boundaries_arr1[t,p,rad] = 2
;                                        if boundaries_arr[t,p,rad] eq 3 then flux[t,p,0]=0.5+flux[t,p,0] ; Trans PUI, kappa0 by MZK
 
					flux_test1[t,p,rad]=calc_diff_flux(dist,energy,r,tp1,ur,ut,nn,np1,theta,kappaPUI,bet[t,p,rad])
					if boundaries_arr[t,p,rad] eq 6 then flux_test2[t,p,rad]=calc_diff_flux(dist,energy,r,tp2,ur,ut,nn,np2,theta,kappaOH,bet[t,p,rad])
					if boundaries_arr[t,p,rad] ne 6 then flux_test2[t,p,rad]=calc_diff_flux(dist,energy,r,tp2,ur,ut,nn,np2,theta,kappaRef,bet[t,p,rad])
					if giacalone ne 1 then flux_test3[t,p,rad]=calc_diff_flux(dist,energy,r,tp3,ur,ut,nn,np3,theta,kappaEn,bet[t,p,rad]) else $
                                                               flux_test3[t,p,rad]=calc_diff_flux_zank(dist,energy,r,ur,ut,nn,np3,theta,kappaEn,bet[t,p,rad],een_frac,em_frac,vshock)
					flux_test_sw[t,p,rad]=calc_diff_flux(dist,energy,r,tsw,ur,ut,nn,nsw,theta,kappa0,bet[t,p,rad])

                        ;Track parameters for globally passed variable
;                                        radial_profile.Temp(0,t,p,rad) = tp & radial_profile.Temp(1,t,p,rad) = tsw &radial_profile.Temp(2,t,p,rad) = tpu3
 ;                                       radial_profile.density(0,t,p,rad) = np & radial_profile.density(1,t,p,rad) = nsw &radial_profile.density(2,t,p,rad) = npu3
                                        radial_profile.nH(t,p,rad) = nn  & radial_profile.Ur(t,p,rad) =data.vr[t,p,rad]^2.  & radial_profile.Ut(t,p,rad) =  data.vt[t,p,rad]
                                        radial_profile.vtheta(t,p,rad) = utheta  & radial_profile.vphi(t,p,rad) = uphi
                                        radial_profile.eram(t,p,rad ) = [.5*mp*data.density[t,p,rad]*(data.vr[t,p,rad]^2. + data.vt[t,p,rad]^2.)*(10^5.)^2.]

                                        radial_profile.emag(t,p,rad) = (ophermagnetic.Bx[t,p,rad]^2. +ophermagnetic.By[t,p,rad]^2. +ophermagnetic.Bz[t,p,rad]^2.)/(8.*!pi)
                        ;------------------- 

; COMMENT OUT LATER - PROXY TEST FOR ACCELERATION IN IHS FOR 8.38 keV
;if eband eq 14 and ext[t,p,rad] gt 1/exp(2) and ext[t,p,rad] lt 1. then begin
;        np1=np1/3.5
;        np2=np2/3.5
;        np3=np3/3.5
;endif

                                endelse
                        endif
                endif

; END OF ADDED STUFF
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
;vplasma=abs(VenaErf+vr)

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
;				vplasma[i,j,k]=ena_frame(Ea,vr[i,j,k])
				vplasma[i,j,k]=plasma_frame2(Ea,vr[i,j,k],vn[i,j,k],tp[i,j,k])
				E_rel[i,j,k]=0.5*mp*(vplasma[i,j,k])^2.
				if E_rel[i,j,k] ne E_rel[i,j,k] then stop
				sigma[i,j,k]=xsection(E_rel[i,j,k])
;                        	sigma[i,j,k]=xsection(Ea,vena)
				bet[i,j,k]=double(sigma[i,j,k]*vplasma[i,j,k]*np[i,j,k]*dr)
;                        	bet[i,j,k]=double(sigma*(vplasma[i,j,k]^2.)*nh[i,j,k]*dr/abs(vr[i,j,k]))
;				if k eq 0 then bet[i,j,k]=bet[i,j,k]+beta0(1.-1./30.) ; considering ionization rate from R=0 to R=30 AU (inner edge of inner boundary)
				bet[i,j,k]=bet[i,j,k]+bet[i,j,k-1]	
			endif
		endfor
	endfor
endfor

end
