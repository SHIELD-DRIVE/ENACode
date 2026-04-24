;=======================================================
;
;     Makes and plots countrate and flux maps 
;
;
;===================================

;Top level program, which plots fluxes for each energy band with mand
;without noise as well as count rates for each modeler and each distribution type 
;-------------------------------------------------------------
pro plot_all

common kappa_value,kappa
common energybands,num_ebands,arr_eband_indexes
common model_names,models,num_models
common distributions,dist_names,num_dist
common plot_flags, onerange_flag, manual_onerange_flag, plot_log, plot_limited, manual_lim, manual_maxset, manual_minset, eps_plot
common plot_params,latmin,latmax,lonmin,lonmax,map_dlat,map_dlon,map_numlatlabels,map_lons,map_lonnames,onerange_flux_min, onerange_flux_max,onerange_cr_min, onerange_cr_max

;most of this code is for calculating over all max and min for each energy band
maximum={plotmax,flux:fltarr(fix(num_ebands)),cr:fltarr(fix(num_ebands))}
minimum={plotmin,flux:fltarr(fix(num_ebands)),cr:fltarr(fix(num_ebands))}
fluxmax=0.
crmax=0.
crmin=0.
fluxmin=0.
minimum.flux[*]=fluxmin
minimum.cr[*]=crmin

for i=0,num_ebands-1 do begin; 
       energyband = arr_eband_indexes(i)
    ;find min and max for each e band for each model and distribution
    if onerange_flag eq 1 then begin
     	for data=0, num_models-1 do begin      
       		for disttype=0, num_dist-1 do begin  
		    read_countrate,energyband,disttype,data,countrate_map
		    make_fluxmap,countrate_map,energyband,flux_map
                    fluxmax = max(flux_map) & fluxmin = min(flux_map)
        	    crmax = max(countrate_map) & crmin = min(countrate_map)
	            if manual_onerange_flag ne 1 then begin
        	   	 if fluxmax gt maximum.flux[i] then maximum.flux[i]=fluxmax
   	           	 if fluxmin lt minimum.flux[i] then minimum.flux[i]=fluxmin	
   	           	 if crmax gt maximum.cr[i] then maximum.cr[i]=crmax
   	            	 if crmin lt minimum.cr[i] then minimum.cr[i]=crmin
		    endif
   		endfor
   	endfor
   endif
   if manual_onerange_flag eq 1 then begin
	maximum.flux[i] = onerange_flux_max[i]
	minimum.flux[i] = onerange_flux_min[i]
	maximum.cr[i] = onerange_cr_max[i]
	minimum.cr[i] = onerange_cr_min[i]
   endif
   for data=0, num_models-1 do begin  
        for disttype=0, num_dist-1 do begin 
   	    read_countrate,energyband,disttype,data,countrate_map
	    make_fluxmap,countrate_map,energyband,flux_map
            fluxmax = max(flux_map) & fluxmin = min(flux_map)
       	    crmax = max(countrate_map) & crmin = min(countrate_map)	            
	    if onerange_flag ne 1 and manual_onerange_flag ne 1 then begin
                  maxflux = fluxmax & minflux = fluxmin
                  maxcr= crmax & mincr = crmin
            endif else begin
            	  maxflux = maximum.flux[i] & minflux = minimum.flux[i]
	          maxcr = maximum.cr[i] & mincr = minimum.cr[i]
            endelse
	    ;make plots
	    plot_averageflux,flux_map,energyband,disttype,data,maxflux,minflux,i 
            ;plot_cr,countrate_map,energyband, disttype,data,maxcr,mincr
        endfor
    endfor
endfor

end
;----------------------------------------------
;reads count rates data files
pro read_countrate,eband,dist_type,model_type,countrate_map

common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher
common kappa_value,kappa
common path_names,countrate_path,averagefluxplots_path,countrateplots_path,count_input_path,count_output_path
common model_names,models,num_models
common distributions,dist_names,num_dist

modeler = models(model_type)
if modeler ne 'Opher_multiion_withneutrals' then countrate_map=fltarr(ntheta,nphi) else countrate_map=fltarr(ntheta,nphi,5)
;if modeler ne 'Opher_multiion_withneutrals' then countrate_map=fltarr(ntheta,nphi) else countrate_map=fltarr(ntheta,nphi,4)

path=countrate_path+count_output_path
dist = dist_names(dist_type)

close,1 ;make sure lun is closed
if modeler ne 'Opher_multiion_withneutrals' then begin
openr,1,path+modeler+'CountRate'+dist+'Eband'+string(eband,format='(i2.2)')+'.dat'
i=0. & j=0.
while not eof(1) do begin
    readf,1,line
    countrate_map[i,j]=line
    j++
    if j eq nphi then begin
        j=0.
        i++
    endif
endwhile
close,1
endif else begin
;for k = 0,3 do begin ; Zank MI
for k=0,4 do begin ; Zank MI + Inj PUIs & Malama MI
;	if k eq 0 then tag = 'swplasma' & if k eq 1 then tag = 'ptplasma' & if k eq 2 then tag = 'prplasma' & if k eq 3 then tag = 'injplasma' & if k eq 4 then tag='combplasma' ; Zank MI + Inj PUIs
;       if k eq 0 then tag = 'swplasma' & if k eq 1 then tag = 'ptplasma' & if k eq 2 then tag = 'prplasma' & if k eq 3 then tag='combplasma' ; Zank MI
       if k eq 0 then tag = 'PUI1' & if k eq 1 then tag = 'PUI2' & if k eq 2 then tag = 'PUI3' & if k eq 3 then tag = 'SW' & if k eq 4 then tag='CombPUIs' 
	openr,1,path+modeler+tag+'CountRate'+dist+'Eband'+string(eband,format='(i2.2)')+'.dat'
i=0. & j=0.
while not eof(1) do begin
    readf,1,line
    countrate_map[i,j,k]=line
    j++
    if j eq nphi then begin
        j=0.
        i++
    endif
endwhile
   close,1
endfor
endelse

end
;--------------------------------------------
pro make_fluxmap,countrate_map,eband,flux_map
;turn countrate back into flux

common sensor_params, hi,lo, wt9, wt10, wt11, wt12, wt13, ultra

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Added by MZK
;
;eminarr=[lo.E.min,hi.E.min]
;emaxarr=[lo.E.max,hi.E.max]
;emin=Eminarr[Eband]
;denergy=(Emaxarr[Eband]-Eminarr[Eband])/10. ; Eres=10
;energy=findgen(10)
;energy=energy*denergy+emin
;avgenergy=total(energy)/10.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ecentral=[lo.e.central,hi.e.central, ultra.e.central]
g=[lo.g,hi.g,ultra.g]
flux_map=countrate_map/ecentral[eband]/g[eband]
print, max(flux_map(*,*,4))
return
end
;--------------------------------------------

pro plot_averageflux,averageflux_map,eband,dist_type,model_type,maxflux,minflux,band_num
;plots average flux over an energy band
;Inputs
;       eband =energyband number 0-7 lo 8-13 hi
;      
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_,rf_opher
common path_names,countrate_path,averagefluxplots_path,countrateplots_path,count_input_path,count_output_path
common plot_flags, onerange_flag, manual_onerange_flag, plot_log, plot_limited, manual_lim, manual_maxset, manual_minset, eps_plot
common plot_params,latmin,latmax,lonmin,lonmax,map_dlat,map_dlon,map_numlatlabels,map_lons,map_lonnames,onerange_flux_min, onerange_flux_max,onerange_cr_min, onerange_cr_max
common model_names,models,num_models
common distributions,dist_names,num_dist
common sensor_params, hi,lo, wt9, wt10, wt11, wt12, wt13, ultra

ecentral=[lo.e.central,hi.e.central,ultra.e.central]
dist = dist_names(dist_type)
path=averagefluxplots_path+count_output_path+'maps/'
modeler = models(model_type)

if modeler ne 'Opher_multiion_withneutrals' then begin
	filename=path+'Eband'+string(eband,format='(i2.2)')+modeler+'averageflux'+dist
if not keyword_set(onerange_flag) and not keyword_set(manual_onerange_flag) then begin
	findmax=max(averageflux_map)
	findmin=min(averageflux_map)
 	maxflux=findmax
 	minflux=findmin
endif

;save fluxes here 
	openw, lun, filename + '.dat', /get_lun
	for i = 0,ntheta-1 do begin
		for j =0,nphi -1 do begin
			printf, lun, i, j, averageflux_map(i,j);ntheta, nphi, value
		endfor
	endfor
	close,lun
	free_lun, lun 

;	title=modeler+', Energy= '+string(ecentral[eband],format='(f5.3)')+' keV, '+dist
	title='Energy = '+string(ecentral[eband],format='(f6.2)')+' keV'
	units='ENA/(cm!U2!Ns sr keV)' 
	if eps_plot eq 1 then plotter_eps,averageflux_map,filename,title,units,maxflux,minflux,band_num else plotter,averageflux_map,filename,title,units,maxflux,minflux,band_num
	
endif else begin
	tempflux_map = fltarr(ntheta,nphi)
;	for k = 0,3 do begin ; Zank MI
	for k=0,4 do begin ; Zank MI + Inj PUIs & Malama MI
;		if k eq 0 then tag = 'swplasma' & if k eq 1 then tag = 'ptplasma' & if k eq 2 then tag = 'prplasma' & & if k eq 3 then tag = 'injplasma' & if k eq 4 then tag = 'combplasma' ; Zank MI + Inj
;               if k eq 0 then tag = 'swplasma' & if k eq 1 then tag = 'ptplasma' & if k eq 2 then tag = 'prplasma' & & if k eq 3 then tag = 'combplasma' ; Zank MI
               if k eq 0 then tag = 'PUI1' & if k eq 1 then tag = 'PUI2' & if k eq 2 then tag = 'PUI3' & & if k eq 3 then tag = 'SW' & if k eq 4 then tag = 'CombPUIs' ; Malama MI
		filename=path+'Eband'+string(eband,format='(i2.2)')+modeler+tag+'averageflux'+dist
	;save fluxes here 
		openw, lun, filename + '.dat', /get_lun
		for i = 0,ntheta-1 do begin
			for j =0,nphi -1 do begin
				printf, lun, i, j, averageflux_map(i,j,k);ntheta, nphi, value
			endfor
		endfor
		close,lun
		free_lun, lun 
;		title=modeler+tag+', Energy= '+string(ecentral[eband],format='(f5.3)')+' keV, '+dist
		title='Energy = '+string(ecentral[eband],format='(f6.2)')+' keV'

		units='ENA(cm!U2!Ns sr keV)' 
		tempflux_map(*,*) = averageflux_map(*,*,k)
		if not keyword_set(onerange_flag) and not keyword_set(manual_onerange_flag) then begin
			findmax=max(tempflux_map)
			findmin=min(tempflux_map)
 			maxflux=findmax
 			minflux=findmin
		endif 
               ; if k eq 2 then stop
		if eps_plot eq 1 then plotter_eps,tempflux_map,filename,title,units,maxflux,minflux,band_num else plotter,tempflux_map,filename,title,units,maxflux,minflux,band_num
	endfor
        ;sum of PUI and SW
        ;tag = 'SW+PU3'
        ;filename=path+'Eband'+string(eband,format='(i2.2)')+modeler+tag+'averageFlux'+dist 
        ;title=modeler+tag+', Energy= '+string(ecentral[eband],format='(f5.3)')+' keV, '+dist
        ;tempflux_map(*,*) = averageflux_map(*,*,1) + averageflux_map(*,*,2)
	;if not keyword_set(onerange_flag) and not keyword_set(manual_onerange_flag) then begin
	;	findmax=max(tempflux_map)	
	;	findmin=min(tempflux_map)
 	;	maxflux=findmax
 	;	minflux=findmin
	;endif
        ;plotter,tempflux_map,filename,title,units,maxflux,minflux,band_num
endelse

;save, averageflux_map, filename = filename + '_fluxmap.sav'

return
end
;----------------------------------------------------------------
pro plot_cr,countrate_map,eband,dist_type,model_type,crmax,crmin
;Plots the countrate data written by the countrate program

common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_,rf_opher
common path_names,countrate_path,averagefluxplots_path,countrateplots_path,count_input_path,count_output_path
common model_names,models,num_models
common distributions,dist_names,num_dist
common sensor_params, hi,lo, wt9, wt10, wt11, wt12, wt13, ultra

ecentral=[lo.e.central,hi.e.central,ultra.e.central]

path=countrateplots_path+count_output_path+'maps/'


dist = dist_names(dist_type)
path=countrateplots_path
modeler = models(model_type)

;filename=path+'Eband'+string(eband,format='(i2.2)')+modeler+'CountRate'+dist
;title=modeler+', '+'Energy='+string(ecentral[eband],format='(f5.3)')+' keV, '+dist
;units='Counts/s'
;plotter,countrate_map,filename,title,units, crmax, crmin

if modeler ne 'Opher_multiion_withneutrals' then begin
	filename=path+'Eband'+string(eband,format='(i2.2)')+modeler+'CountRate'+dist
;	title=modeler+', Energy= '+string(ecentral[eband],format='(f5.3)')+' keV, '+dist
	title='Energy = '+string(ecentral[eband],format='(f6.3)')+' keV'
	units='Counts/s'
	if not keyword_set(onerange_flag) and not keyword_set(onerange_manual_flag) then begin
			findmax=max(countrate_map)
			findmin=min(countrate_map)
 			maxflux=findmax
 			minflux=findmin
	endif
	plotter,countrate_map,filename,title,units,maxflux,minflux
endif else begin
	tempcountrate_map = fltarr(ntheta,nphi)
	for k = 0,2 do begin
		if k eq 0 then tag = 'totplasma' & if k eq 1 then tag = 'swplasma' & if k eq 2 then tag = 'pu3plasma'
		filename=path+'Eband'+string(eband,format='(i2.2)')+modeler+tag+'CountRate'+dist 
;		title=modeler+tag+', Energy= '+string(ecentral[eband],format='(f5.3)')+' keV, '+dist
		title='Energy = '+string(ecentral[eband],format='(f6.3)')+' keV'
		units='Counts/s' 
		tempcountrate_map(*,*) = countrate_map(*,*,k)
		if not keyword_set(onerange_flag) and not keyword_set(onerange_manual_flag) then begin
			findmax=max(tempcountrate_map)
			findmin=min(tempcountrate_map)
 			maxflux=findmax
 			minflux=findmin
		endif
		plotter,tempcountrate_map,filename,title,units,maxflux,minflux
	endfor
        ;sum of PUI and SW
	tag = 'SW+PU3'
        filename=path+'Eband = '+string(eband,format='(i2.2)')+modeler+tag+'CountRate'+dist 
;        title=modeler+tag+', Energy= '+string(ecentral[eband],format='(f5.3)')+' keV, '+dist
	title='Energy'+string(ecentral[eband],format='(f6.3)')+' keV'
        tempcountrate_map(*,*) = countrate_map(*,*,1) + countrate_map(*,*,2)
	if not keyword_set(onerange_flag) and not keyword_set(onerange_manual_flag) then begin
		findmax=max(tempcountrate_map)	
		findmin=min(tempcountrate_map)
 		maxflux=findmax
 		minflux=findmin
	endif
        plotter,tempcountrate_map,filename,title,units,maxflux,minflux,band_num
endelse




end

