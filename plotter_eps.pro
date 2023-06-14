pro plotter_eps, map,filename, map_title, units, maxset, minset, band_num
common plot_params,latmin,latmax,lonmin,lonmax,map_dlat,map_dlon,map_numlatlabels,map_lons,map_lonnames, onerange_flux_min, onerange_flux_max, onerange_cr_min, onerange_cr_max
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen
common plot_flags, onerange_flag,manual_onerange_flag, plot_log,plot_limited, manual_lim, manual_maxset, manual_minset, eps_plot, plot_scl,tail,spec_min, spec_max, pres_min, pres_max, write_ecl,legend_style
common model_names,models,num_models,pui_model,vasyliunas,ribbon, gdf, hinterp,ntr_frac,nref_frac,ttr_frac,tref_frac,colde,esw_frac,etr_frac,eref_frac,moscow

; flips array upside down (so that it plots 0 theta at the top of the map)
map1=fltarr(nphi,ntheta)  ;flipped for mappings nomenclature

print, "ENERGY BAND: ", map_title

if gdf eq 0 then begin

;print, filename, map[7,0]*plot_scl
;print, filename, map[5,3]*plot_scl
;print, filename, map[21,2]*plot_scl
print, "Nose flux", map[14,30]
print, "Tail flux", map[14,0]
print, "Flank flux", map[14,45]
print, "Lobe flux", map[7,0]
print, "V2 flux", map[19,29]
print, "Max Flux", max(map*plot_scl)
print, "Min Flux", min(map*plot_scl)

;map[5:6,42:46]=100. ribbon kink in 4.3 keV

map_test=fltarr(2,15,15)
map_test[0,0:14,0:14]=map[7:21,23:37]
map_test[1,0:14,0:6]=map[7:21,53:59]
map_test[1,0:14,7:14]=map[7:21,0:7]

print, "Tail/Nose", mean(map_test[1,*,*]/map_test[0,*,*])

; lobe test
;map[24:28,56:59]=0.
;map[24:28,0:5]=0.

; N. Lobe
;map[3:10,57]=0.
;map[3,57:59]=0.
;map[3,0:9]=0.
;map[3:10,9]=0.
;map[10,57:59]=0.
;map[10,0:9]=0.

; Lobe Centers
;map[26,0]=0.
;map[4,2]=0.


center=0
if center eq 1 then begin
   ti=24
   pi=56
   tf=28
   pf=5
   dt=tf-ti
   if pf lt pi then dp=(pf+60)-pi else dp=pf-pi
   tot=0.
   tot_i=0.
   tot_j=0.
   for i=0,dp do for j=0,dt do tot=tot+map[ti+j,pi-60+i]
   for i=0,dp do for j=0,dt do tot_i=map[ti+j,pi-60+i]*i+tot_i
   for i=0,dp do for j=0,dt do tot_j=map[ti+j,pi-60+i]*j+tot_j
   tc=ti+round(tot_j/tot)
   pc=pi+round(tot_i/tot)
   if pc gt 59 then pc=pc-60
   print, tc, pc  
   map[tc,pc]=0.
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; 

;if band_num eq 0 then map=map/0.242529
;if band_num eq 1 then map=map/0.0392021
;if band_num eq 2 then map=map/0.0120232
;if band_num eq 3 then map=map/0.00329506

map1tmp=map*plot_scl
;map1tmp=map/min(map)
;map1tmp=map/max(map)

for i=0,ntheta-1 do begin
    for j=0, nphi-1 do begin
        map1[j,(ntheta-1-i)]=map1tmp[i,j]
    end
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; map2tmp removed since no longer needed
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Reorder grid
map3tmp = map1
for j=0, nphi-1 do begin
        map1[nphi-1-j,*]=map3tmp[j,*]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Making tail at the center of the maps - MZK Dec. 11 2015

nflip=nphi/2

; For having downwind at center when using outputs from ENA Code
map4atmp=map1
if write_ecl eq 0 then begin
   for j=0, nflip-1 do begin
      map1[j,*]=map4atmp[nflip+j,*]
      map1[nflip+j,*]=map4atmp[j,*]
   endfor
endif

map4btmp=map1                    ; changing outside-in to inside-out for eps
for j=nphi-1, 0, -1 do map1[nphi-1-j,*]=map4btmp[j,*]

;map4ctmp=map1

;for j=ntheta-1, 0, -1 do map1[*,ntheta-1-j]=map4ctmp[*,j]; latitudinal flip as proxy for BISM flip

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Convert from simulation to J2000 coords - MZK Nov. 2019

latcon=round(5./dtheta) ; - latitude
loncon=round(281./dphi) ; + longitude

map5tmp = map1 ; For longitude change

for j=0, loncon-1 do begin
   if (loncon+j) le nphi-1 then map1[j,*]=map5tmp[loncon+j,*]
   if (loncon+j) gt nphi-1 then map1[j,*]=map5tmp[j-round(nphi-loncon),*]
endfor

for j=0, nphi-loncon-1 do map1[loncon+j,*]=map5tmp[j+loncon-round(nphi-loncon),*]

map6tmp = map1 ; For latitude change

;for j=0, ntheta-latcon-1 do map1[*,j]=map6tmp[*,j+latcon]
;map1[*,ntheta-latcon]=map6tmp[*,ntheta-latcon] ; Saying n. pole equals 6 degrees based on nature of how we interpret poles in code

; Flip map latitudes
;for i=0,ntheta-1 do map1[*,i]=map6tmp[*,ntheta-1-i]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if write_ecl ne 0 then begin
tag='Ecl'
openw,1,'MoscowEclipticEband'+string(band_num,format='(i2.2)')+'.dat'
  for j =0, nphi-1 do begin
      for i=0, ntheta-1 do begin
          printf,1,3.+j*6., -87.+i*6., map1[j,i]/plot_scl
      endfor
  endfor
close,1
endif

endif

if ribbon eq 1 then begin
	; Read in ribbon files
	map_1keV=make_array(30,60)
	map_2_7keV=make_array(30,60)
	map_4keV=make_array(30,60)

	file_1keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/f2-ribbonMaps/hi10-rib.txt'
	file_2_7keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/f2-ribbonMaps/hi12-rib.txt'
	file_4keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/f2-ribbonMaps/hi13-rib.txt'

	map_1keV=read_table(file_1keV)
	map_2_7keV=read_table(file_2_7keV)
        map_4keV=read_table(file_4keV)

        map_1keV_tmp=map_1keV
        map_2_7keV_tmp=map_2_7keV
        map_4keV_tmp=map_4keV

        for j=0, nflip-1 do begin
           map_1keV[j,*]=map_1keV_tmp[nflip+j,*]
           map_2_7keV[j,*]=map_2_7keV_tmp[nflip+j,*]
           map_4keV[j,*]=map_4keV_tmp[nflip+j,*]

           map_1keV[nflip+j,*]=map_1keV_tmp[j,*]
           map_2_7keV[nflip+j,*]=map_2_7keV_tmp[j,*]
           map_4keV[nflip+j,*]=map_4keV_tmp[j,*]
        endfor
           
	; Add ribbon flux to map flux
	if band_num eq 0 then map1=map1*5+map_1keV
	if band_num eq 1 then map1=map1*15.+map_2_7keV
	if band_num eq 2 then map1=map1*60.+map_4keV
endif

if moscow eq 1 then begin

;   	moscow_file=read_ascii('MoscowMHD_ZirnsteinModel_ENAdata.dat', data_start=1)
;        moscow_file=read_ascii('alpha_PUI_TS.dat', data_start=1)
        moscow_file=read_ascii('MoscowEclipticEband00.dat')
   	map_0_7keV=make_array(60,30)
   	map_1keV=make_array(60,30)
        map_1_7keV=make_array(60,30)
        map_2_7keV=make_array(60,30)
        map_4keV=make_array(60,30)

        for j=0, 59 do begin 
		for i=29,0,-1 do begin
			map_0_7keV(j,i)=moscow_file.field1(2,i+j*30)
                        map_1keV(j,i)=moscow_file.field1(3,i+j*30)
                        map_1_7keV(j,i)=moscow_file.field1(4,i+j*30)
                        map_2_7keV(j,i)=moscow_file.field1(5,i+j*30)
                        map_4keV(j,i)=moscow_file.field1(6,i+j*30)
		endfor
	endfor

        map_0_7keV_tmp=map_0_7keV
        map_1keV_tmp=map_1keV
        map_1_7keV_tmp=map_1_7keV
        map_2_7keV_tmp=map_2_7keV
        map_4keV_tmp=map_4keV

        nflip=30.; for 6 degree res. of IBEX, have 60 long. bins (nphi/2)

        for j=0, nflip-1 do begin
           map_0_7keV[j,*]=map_0_7keV_tmp[nflip+j,*]
           map_0_7keV[nflip+j,*]=map_0_7keV_tmp[j,*]
           map_1keV[j,*]=map_1keV_tmp[nflip+j,*]
           map_1keV[nflip+j,*]=map_1keV_tmp[j,*]
           map_1_7keV[j,*]=map_1_7keV_tmp[nflip+j,*]
           map_1_7keV[nflip+j,*]=map_1_7keV_tmp[j,*]
           map_2_7keV[j,*]=map_2_7keV_tmp[nflip+j,*]
           map_2_7keV[nflip+j,*]=map_2_7keV_tmp[j,*]
           map_4keV[j,*]=map_4keV_tmp[nflip+j,*]
           map_4keV[nflip+j,*]=map_4keV_tmp[j,*]
        endfor

        ; Add ribbon flux to map flux
        if band_num eq 0 then map1=map_0_7keV*plot_scl
        if band_num eq 1 then map1=map_1keV*plot_scl
        if band_num eq 2 then map1=map_1_7keV*plot_scl
        if band_num eq 3 then map1=map_2_7keV*plot_scl
        if band_num eq 4 then map1=map_4keV*plot_scl


endif

if gdf eq 1 then begin

        ; IBEX
        ; Read in GDF files
        map_0_7keV=make_array(30,60)
        map_1keV=make_array(30,60)
   	map_1_7keV=make_array(30,60)
	map_2_7keV=make_array(30,60)
	map_4keV=make_array(30,60)

	sn_0_7keV=make_array(30,60)
        sn_1keV=make_array(30,60)
        sn_1_7keV=make_array(30,60)
        sn_2_7keV=make_array(30,60)
        sn_4keV=make_array(30,60)

	; Schwadron+2014 GDF (2009-2013)
	;file_0_7keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/f3-gdfMaps/hi09-gdf.txt'
        ;file_1keV='/home/kmarc/Downloads/RibbonSeparationPaperReleaseDataDirectories/f3-gdfMaps/hi3-gdf.txt'
        ;file_1_7keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/f3-gdfMaps/hi11-gdf.txt'
        ;file_2_7keV='/home/kmarc/Downloads/RibbonSeparationPaperReleaseDataDirectories/f3-gdfMaps/hi5-gdf.txt'
        ;file_4keV='/home/kmarc/Downloads/RibbonSeparationPaperReleaseDataDirectories/f3-gdfMaps/hi6_gdf.txt'

	;file_sn_0_7keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/uncert_gdfMaps/hi09-sn.txt'
	;file_sn_1keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/uncert_gdfMaps/hi10-sn.txt'
	;file_sn_1_7keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/uncert_gdfMaps/hi11-sn.txt'
	;file_sn_2_7keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/uncert_gdfMaps/hi12-sn.txt'
	;file_sn_4keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/uncert_gdfMaps/hi13-sn.txt'

	; McComas+2020 IBEX Data Release 16 (2009-2012)
        file_0_7keV='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/IBEX_2009_2012/ibex_0.71-flux.txt'
        file_1keV='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/IBEX_2009_2012/ibex_1.11-flux.txt'
        file_1_7keV='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/IBEX_2009_2012/ibex_1.74-flux.txt'
        file_2_7keV='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/IBEX_2009_2012/ibex_2.73-flux.txt'
        file_4keV='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/IBEX_2009_2012/ibex_4.29-flux.txt'

        file_sn_0_7keV='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/IBEX_2009_2012/ibex_0.71-fsnr.txt'
        file_sn_1keV='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/IBEX_2009_2012/ibex_1.11-fsnr.txt'
        file_sn_1_7keV='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/IBEX_2009_2012/ibex_1.74-fsnr.txt'
        file_sn_2_7keV='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/IBEX_2009_2012/ibex_2.73-fsnr.txt'
        file_sn_4keV='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/IBEX_2009_2012/ibex_4.29-fsnr.txt'

        map_0_7keV=read_table(file_0_7keV)
	map_1keV=read_table(file_1keV)
        map_1_7keV=read_table(file_1_7keV)
	map_2_7keV=read_table(file_2_7keV)
        map_4keV=read_table(file_4keV)

	sn_0_7keV=read_table(file_sn_0_7keV)
        sn_1keV=read_table(file_sn_1keV)
        sn_1_7keV=read_table(file_sn_1_7keV)
        sn_2_7keV=read_table(file_sn_2_7keV)
        sn_4keV=read_table(file_sn_4keV)

        map_0_7keV_tmp=map_0_7keV
        map_1keV_tmp=map_1keV
        map_1_7keV_tmp=map_1_7keV
        map_2_7keV_tmp=map_2_7keV
        map_4keV_tmp=map_4keV

        sn_0_7keV_tmp=sn_0_7keV
        sn_1keV_tmp=sn_1keV
        sn_1_7keV_tmp=sn_1_7keV
        sn_2_7keV_tmp=sn_2_7keV
        sn_4keV_tmp=sn_4keV

        nflip=30.; for 6 degree res. of IBEX, have 60 long. bins (nphi/2)

        for j=0, nflip-1 do begin
           map_0_7keV[j,*]=map_0_7keV_tmp[nflip+j,*]
           map_0_7keV[nflip+j,*]=map_0_7keV_tmp[j,*]
           map_1keV[j,*]=map_1keV_tmp[nflip+j,*]
           map_1keV[nflip+j,*]=map_1keV_tmp[j,*]
           map_1_7keV[j,*]=map_1_7keV_tmp[nflip+j,*]
           map_1_7keV[nflip+j,*]=map_1_7keV_tmp[j,*]
           map_2_7keV[j,*]=map_2_7keV_tmp[nflip+j,*]
           map_2_7keV[nflip+j,*]=map_2_7keV_tmp[j,*]
           map_4keV[j,*]=map_4keV_tmp[nflip+j,*]
           map_4keV[nflip+j,*]=map_4keV_tmp[j,*]
        endfor
        
        for j=0, nflip-1 do begin
           sn_0_7keV[j,*]=sn_0_7keV_tmp[nflip+j,*]
           sn_0_7keV[nflip+j,*]=sn_0_7keV_tmp[j,*]
           sn_1keV[j,*]=sn_1keV_tmp[nflip+j,*]
           sn_1keV[nflip+j,*]=sn_1keV_tmp[j,*]
           sn_1_7keV[j,*]=sn_1_7keV_tmp[nflip+j,*]
           sn_1_7keV[nflip+j,*]=sn_1_7keV_tmp[j,*]
           sn_2_7keV[j,*]=sn_2_7keV_tmp[nflip+j,*]
           sn_2_7keV[nflip+j,*]=sn_2_7keV_tmp[j,*]
           sn_4keV[j,*]=sn_4keV_tmp[nflip+j,*]
           sn_4keV[nflip+j,*]=sn_4keV_tmp[j,*]
        endfor

	; Add ribbon flux to map flux
        if band_num eq 0 then map1=map_0_7keV
	if band_num eq 1 then map1=map_1keV
        if band_num eq 2 then map1=map_1_7keV
	if band_num eq 3 then map1=map_2_7keV
	if band_num eq 4 then map1=map_4keV

endif

if gdf eq 2 then begin

        ; INCA
        ; Read in GDF files
        file_14='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/INCAData/TOF72009_2012.sav'
        file_15='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/INCAData/TOF62009_2012.sav'
        file_16='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/INCAData/TOF52009_2012.sav'
        file_17='/nfs/magnetic6/kmarc/ENAMAPCODE_KMHD/INCAData/TOF42009_2012.sav'
      
        if band_num eq 0 then restore, file_14
        if band_num eq 1 then restore, file_15
        if band_num eq 2 then restore, file_16
        if band_num eq 3 then restore, file_17
        
        map_tmp=map
        longitude=longitude
        latitude=latitude
        
        nflip=160.; for 6 degree res. of IBEX, have 60 long. bins (nphi/2)

        ;map1=fltarr(nphi,ntheta)
        for j=0, nflip-1 do begin
           map1[j,*]=map_tmp[nflip+j,*]
           map1[nflip+j,*]=map_tmp[j,*]
        endfor

        wmax=where(map1 gt manual_maxset(band_num))
        map1[wmax]=0.    
 
        ;map1[0:floor(nphi/2.),*]=0.
        ;map1[nphi-floor(nphi/2.):nphi-1,*]=0.
 
endif

; S. Lobe @ [45,8] when using 6 deg x 6 deg resolution
; V1 @ [17,21] when using 6 deg x 6 deg resolution
; V2 @ [12,10] when using 6 deg x 6 deg resolution
; Using 18 deg x 18 deg blocks centered on location
if gdf eq 1 then begin 
		print, filename, "  S. Lobe", mean(map1[44:46,7:9]), " +/-", mean(map1[44:46,7:9]/sn_0_7keV[44:46,7:9])
        	print, filename, "  V1", mean(map1[11:13,20:22]), " +/-", mean(map1[11:13,20:22]/sn_0_7keV[11:13,20:22])
                print, filename, "  V2", mean(map1[16:18,9:11]), " +/-", mean(map1[16:18,9:11]/sn_4keV[16:18,9:11])
                print, filename, "  Flank", mean(map1[27:29,14:16]), " +/-", mean(map1[27:29,14:16]/sn_4kev[27:29,14:16])
                print, filename, "  Flank 2", mean(map1[23:25,14:16]), " +/-", mean(map1[23:25,14:16]/sn_4kev[23:25,14:16])
                print, filename, "  Tail", mean(map1[42:44,14:16]), " +/-", mean(map1[42:44,14:16]/sn_0_7kev[42:44,14:16])
                print, filename, "  Nose", mean(map1[12:14,13:15]), " +/-", mean(map1[42:44,14:16]/sn_0_7kev[12:14,13:15])
endif else if gdf eq 0 then begin
	print, filename, "  S. Lobe",  mean(map1[44:46,7:9])
        print, filename, "  V1", mean(map1[11:13,20:22])
        print, filename, "  V2", mean(map1[16:18,9:11])
        print, filename, "  Flank", mean(map1[27:29,14:16]) ; 90 deg from Nose (like Giacalone+2021)
        print, filename, "  Flank 2", mean(map1[23:25,14:16]) ; Avoids NaN snr in IBEX 2009-2012 data
        print, filename, "  Tail", mean(map1[42:44,14:16])
        print, filename, "  Nose", mean(map1[12:14,13:15])

;print, filename, "  Between Lobe", mean(map1[65:67,22:23])
;print, filename, "  N. Lobe Peak", mean(map1[69:71,26:27])
;print, filename, "  N. Lobe Peak/Between Lobe", mean(map1[68:70,29:30])/mean(map1[65:67,22:23])
endif

;map1[*,*]=9999.
;map1[16:18,9:11]=0.
;map1[23:25,14:16]=0.
;map1[42:44,14:16]=0.
;map1[65:67,22:23]=0. ; 4 deg between lobes
;map1[69:71,26:27]=0. ; 4 deg moscow n. lobe peak (80 kev)
;map1[68:70,29:30]=0. ; 4 deg moscow n. lobe peak (44 kev)
;map1[70:72,24:25]=0. ; 4 deg bu n. lobe peak (80 kev)
;map1[69:71,28:29]=0. ; 4 deg bu n. lobe peak (44 kev)

if plot_limited eq 1 then begin
        plot_limits =[latmin,lonmin,latmax,lonmax]
        ;plot_limits =[0,180,90,0,0,-180,-90,0]
        ;if not in range, don't count as max
        x_lo = floor((nphi-1)/2.+fix(lonmin/dphi)) & x_hi = floor((nphi-1)/2.+fix(lonmax/dphi))
        y_lo = floor((ntheta-1)/2.+fix(latmin/dtheta)) & y_hi = floor((ntheta-1)/2.+fix(latmax/dtheta))
	if x_lo eq -1 then x_lo = 0 & if y_lo eq -1 then y_lo = 0; special case 
	if onerange_flag ne 1 and manual_onerange_flag ne 1 then maxset = max(map1[x_lo:x_hi,y_lo:y_hi])
endif


;makes sure all of the values are positive
for i=0, nphi-1 do begin
    for j=0, ntheta-1 do begin
        if map1[i,j] lt 0.0 then map1[i,j]=0.0
    endfor
endfor

;scale colors
if plot_log ne 1 then begin
    ;changes map into a color map
    if manual_lim eq 1 then begin ; if using manual color bar limits -MZK (24 May 2016) 
	maxset=manual_maxset(band_num) ; manual max flux value for maps
	minset=manual_minset(band_num) ; manual min flux value for maps
	for i=0,nphi-1 do begin
		for j=0,ntheta-1 do begin
			if map1[i,j] lt minset then map1[i,j]=minset
		endfor
	endfor
    endif
    map2=254.*(map1-minset)/(maxset-minset)
    above = where(map2 gt 254.)
    above_size = size(above)
    if above_size(0) ne 0 then map2(above) = 253.
    if manual_lim eq 1 then maxcolor=max(map2) ; if using manual plot limits - MZK
    if manual_lim ne 1 then maxcolor = 253 ; if using default plot limits - MZK
endif

if plot_log eq 1 then  begin
    if manual_lim eq 1 then begin ; if using manual color bar limits -MZK (24 May 2016) 
        maxset=manual_maxset(band_num) ; manual max flux value for maps
        minset=manual_minset(band_num) ; manual min flux value for maps
	for i=0,nphi-1 do begin
                for j=0,ntheta-1 do begin
                        if map1[i,j] lt minset then map1[i,j]=minset
                endfor
        endfor
    endif
    map2=254.*(alognum(map1)-alognum(minset))/(alognum(maxset)-alognum(minset))
    ;map2 = 255.*((map1-minset)/(maxset-minset))^(1./2.)
    ;sqrt scheme
    if manual_lim eq 1 then maxcolor=max(map2) ; if using manual plot limits - MZK
    if manual_lim ne 1 then maxcolor = 253 ; if using default plot limits - MZK
endif

window,0,retain=2 ; gets plot formatting correct without jpeg

xsize=9.0
ysize=6.0
dpi=100
margin=0.0
ps_file=filename+'.eps'
ps_on, filename=ps_file, margin=margin, page_size=[xsize,ysize], /inches
device,decomposed=0
!p.background = 255.
!p.color = 0.
cgloadct,33,/silent,ncolors=254
;cgloadct,70,/silent,ncolors=254
;cgloadct,0,/silent,ncolors=254
if tail eq 1 then begin
	if plot_limited eq 1 then cgmap_set,-5,79,/HAMMER,/ISOTROPIC,/noerase, color=0, charsize = 2.5,position=[0.0148, 0.214, 0.983, 0.942],limit = plot_limits, /noborder else map_set,-5,79,/HAMMER,/ISOTROPIC, color=0,charsize = 2.5, position=[0.0148, 0.214, 0.983, 0.942],/noborder, reverse=1
endif else begin
	if plot_limited eq 1 then cgmap_set,5,-101,/HAMMER,/ISOTROPIC,/noerase, color=0, charsize = 2.5,position=[0.0148, 0.214, 0.983, 0.942],limit = plot_limits, /noborder else map_set,5,-101,/HAMMER,/ISOTROPIC, color=0,charsize = 2.5, position=[0.0148, 0.214, 0.983, 0.942],/noborder, reverse=1
endelse
if gdf ge 1 then result=map_image(map2,startx,starty,missing=255,compress=1) else $
result = map_image(map2,startx,starty,missing=255,compress=1);,/bilinear)
thisposition=[0.02, 0.212, 0.97, 0.94] ;[0.0148, 0.212, 0.983, 0.94];[0.0148,0.112,0.983,0.84]
cgimage,result,position=thisposition

; Overplotting BV plane onto ENA Maps
;map3=make_array(360,180)
;helio2ecl, map3
;;mapbvtmp=map3
;;nt=180.
;;for j=nt-1, 0, -1 do map3[*,nt-1-j]=mapbvtmp[*,j]; latitudinal flip as proxy for BISM flip
;result2 = map_image(map3,startx,starty,missing=0,/bilinear,compress=1)
;nwidth=n_elements(result2[*,0])
;nheight=n_elements(result2[0,*])
;for t=0,1 do begin
;   for i=0,nwidth-1 do begin
;      for j=0,nheight-1 do begin
;         if t eq 0 then if result2[i,j] ne 0 then result2[i,j]=255.
;         if t eq 1 then if result2[i,j] eq 0 then result2[i,j]=result[i,j]
;      endfor
;   endfor
;endfor
;cgimage,result2,position=thisposition

cgloadct,39,/silent
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plotting Energy on Map
if plot_scl eq 1 then xyouts,350,14505,map_title,/device,color=0,charsize=4.0,charthick=7.
if plot_scl ne 1 then xyouts,350,14505,map_title+' (Flux x '+strtrim(plot_scl,2)+')',/device,color=0,charsize=4.0,charthick=7.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SCALE:
x0 = 70 & y0 = 40 ; starting pos of bar [pix]
xlen = 500 & ylen = 20 ; len of bar [pix]
val = findgen(xlen)/(xlen-1.0)*255.0
bar = fltarr(xlen,ylen)
for i=0,ylen-1 do begin
    bar[*,i] = val
end

STICKN = 5

lats=[-90,-60,-30,0,30,60,90]
lons=[-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180]
map_grid,color=0,glinethick=5,glinestyle=0,label=1,LATS=lats,LONS=lons,charsize=3,charthick=8 ; Labeled grid lines
if tail eq 1 then lons2=[-105,-95] else lons2=[76,84] ; for edges
map_grid,color=0,glinethick=5,glinestyle=0,LONS=lons2,charsize=3,charthick=8 ; Grid lines for edges

;plots voyager and nose positions, correct for nose rotated to 0 (should be 255)
noselong=-101 ;0
taillong=79.;-(75.-255.) ;180
voyager1long = -106;-7.7
voyager2long=-71;34
noselat = 5 ;just above grid lines
taillat= -5.;0
voyager1lat= 35
voyager2lat=-32

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Taken out by MZK on Dec. 11 2015
;
if tail eq 1 then begin
        x= [taillong];
        y=[taillat];
        label=["DOWNWIND"];
endif else begin
        x= [noselong, voyager1long,voyager2long]
        y=[noselat, voyager1lat,voyager2lat]
        label=["UPWIND", "V1","V2"]
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Taken out by MZK on Dec. 11 2015
;
;put specific points with labels on map
if tail eq 1 then begin
    	oplot, [x-.1,x,x+.1], [y-.1,y,y+.1],symsize=2, psym=8, color=0 ; dot for downwind location
    	xyouts,x-3.0,y-10, label, color=0,charsize=4, charthick=10 ; label of downwind
endif else begin
        for i=0, 2 do begin
           oplot, [x[i]-.1,x[i],x[i]+.1], [y[i]-.1,y[i],y[i]+.1],symsize=2, psym=8, color=0
           xyouts,x[i]+4.0,y[i]+5.0, label[i], color=0,charsize=4,charthick=10
        endfor
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
!X.TICKFORMAT = ''
cgloadct,33,/silent,ncolors=254
;cgloadct,70,/silent,ncolors=254

;make colorbar
if plot_log ne 1 then begin
	if legend_style eq 0 then if units eq 'Counts/s' then cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(f10.3)', position =[0.69,0.93,0.99,0.99], charthick=4, charsize =3 else cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(I4)', position =[0.1,0.17,0.9,0.2], charthick=8, charsize =3. ;( colorbar.pro at dfanning.com)
       if legend_style eq 1 then if units eq 'Counts/s' then cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(f10.3)', position =[0.69,0.93,0.99,0.99], charthick=4, charsize =3 else cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(f10.2)', position =[0.1,0.17,0.9,0.2], charthick=8, charsize =3. ;( colorbar.pro at dfanning.com)
       if legend_style eq 2 then if units eq 'Counts/s' then cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(f10.3)', position =[0.69,0.93,0.99,0.99], charthick=4, charsize =3 else cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(e9.1)', position =[0.1,0.17,0.9,0.2], charthick=8, charsize =3. ;( colorbar.pro at dfanning.com)
endif else begin
    ;cb_ticknames =string((maxset/minset)^(findgen(5)/4.)*minset,format='(I3)')
    cb_ticknames =string((maxset-minset)*(findgen(5)/4.)^2.+minset,format='(I3)')

	if units eq 'Counts/s' then cgcolorbar, bottom = 0, divisions = 4., minor =0., ticknames= cb_ticknames, format = '(f10.3)', position =[0.69,0.93,0.99,0.99], charthick=2, charsize =2 else cgcolorbar, bottom = 0, divisions = 4., minor =0., ticknames= cb_ticknames, format = '(I4)', position =[0.75,0.945,0.99,0.995], charthick=3, charsize =2.5  ;( colorbar.pro at dfanning.com)

endelse

;CONTOUR, result,c_thick=2,xticklayout=1,yticklayout=1,xtickname=REPLICATE(' ', 10),ytickname=REPLICATE(' ', 10),POSITION=thisPosition, $
;/NOERASE, XSTYLE=1,YSTYLE=1, yrange=[0,221], xrange=[0,442],LEVELS=[80,143];,LEVELS=[167,185];, NLEVELS=10

ps_off

;save
;image = TVRD(0,0,!D.x_size,!D.y_size,true=3)
;WRITE_JPEG, filename+'.jpg', image, QUALITY=100, TRUE=3

print, ""

end
;-------------------------------------------
FUNCTION ALOGNUM, x
common log_def, log_num
      RETURN, ALOG(x) / ALOG(log_num)
END

