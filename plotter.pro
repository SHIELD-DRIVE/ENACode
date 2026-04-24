pro plotter, map,filename, map_title, units, maxset, minset, band_num
common plot_params,latmin,latmax,lonmin,lonmax,map_dlat,map_dlon,map_numlatlabels,map_lons,map_lonnames, onerange_flux_min, onerange_flux_max, onerange_cr_min, onerange_cr_max
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen
common plot_flags, onerange_flag,manual_onerange_flag, plot_log,plot_limited, manual_lim, manual_maxset, manual_minset,eps_plot, plot_scl,tail
common model_names,models,num_models, colde

; flips array upside down (so that it plots 0 theta at the top of the map)
map1=fltarr(nphi,ntheta)  ;flipped for mappings nomenclature

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; 
map1tmp=map*plot_scl

for i=0,ntheta-1 do begin
    for j=0, nphi-1 do begin
        map1[j,(ntheta-1-i)]=map1tmp[i,j]
    end
end

; Added by MZK to flip about x-axis; comment out above for use 
;map1tmp=map
;for i=0,ntheta-1 do begin
;    for j=0, nphi-1 do begin
;        map1[j,i]=map1tmp[i,j]
;    end
;end
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;THIS REQUIRES SPECIAL CARE BECAUSE ROTATING NOSE to 0 degrees for plotting
;Nose starts at 75+.5*dphi (how cells are split) degrees in map phi, calc which index that is and rotate
;WHY 78?! 255-180 = 75 . . .  plus cell shift, just order in which put into map
;magicnumber = fix((75.+.5*dphi)/360.*nphi)
;map2tmp=map1
;for i=0,ntheta-1 do begin
;    for j=0, nphi-1 do begin
;       if j ge magicnumber then map1[j-(magicnumber),i]=map2tmp[j,i]
;;       if j lt magicnumber then map1[j+(nphi-magicnumber),i]=map2tmp[j,i]
 ;   end
;end

;flip longitude (x axis) so that larger values are at the beginning of the array (for sky map)
;NEED TO MANUALLY fill in longitude lines due to this calculation

map3tmp = map1
for j=0, nphi-1 do begin
        map1[nphi-1-j,*]=map3tmp[j,*]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Making tail at the center of the maps - MZK Dec. 11 2015

; For having downwind at center when using outputs from ENA Code
if tail eq 1 then begin
	map4tmp=map1
	nflip=nphi/2
	for j=0, nflip-1 do begin
       		map1[j,*]=map4tmp[nflip+j,*]
       		map1[nflip+j,*]=map4tmp[j,*]
	endfor
endif

;map4ctmp=map1

;for j=ntheta-1, 0, -1 do map1[*,ntheta-1-j]=map4ctmp[*,j]; latitudinal flip as proxy for BISM flip

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if ribbon eq 1 then begin
	; Read in ribbon files
	map_1keV=make_array(30,60)
	map_2_7keV=make_array(30,60)
	map_4keV=make_array(30,60)

	;file_1keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/f2-ribbonMaps/hi10-rib.txt'
	;file_2_7keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/f2-ribbonMaps/hi12-rib.txt'
        ;file_4keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/f2-ribbonMaps/hi13-rib.txt'

        file_1keV='/home/kmarc/Downloads/RibbonSeparationPaperReleaseDataDirectories/f3-gdfMaps/hi3-gdf.txt'
        file_2_7keV='/home/kmarc/Downloads/RibbonSeparationPaperReleaseDataDirectories/f3-gdfMaps/hi5-gdf.txt'
        file_4keV='/home/kmarc/Downloads/RibbonSeparationPaperReleaseDataDirectories/f3-gdfMaps/hi6_gdf.txt'

	map_1keV=read_table(file_1keV)
	map_2_7keV=read_table(file_2_7keV)
        map_4keV=read_table(file_4keV)

        map_1keV_tmp=map_1keV
        map_2_7keV_tmp=map_2_7keV
        map_4keV_tmp=map_4keV

        
        ;for j=0, nflip-1 do begin
        ;   map_1keV[j,*]=map_1keV_tmp[nflip+j,*]
        ;   map_2_7keV[j,*]=map_2_7keV_tmp[nflip+j,*]
        ;   map_4keV[j,*]=map_4keV_tmp[nflip+j,*]

        ;   map_1keV[nflip+j,*]=map_1keV_tmp[j,*]
        ;   map_2_7keV[nflip+j,*]=map_2_7keV_tmp[j,*]
        ;   map_4keV[nflip+j,*]=map_4keV_tmp[j,*]
        ;endfor

        nflip=45.-1.
        map_1keV_tmp=make_array(60,30)
        map_1keV_tmp[0:15,*]=map_1keV[nflip:nflip+15,*]
        map_1keV_tmp[16:15+nflip,*]=map_1keV[0:nflip-1,*]
        map_1keV[59:0:-1,*]=map_1keV_tmp[0:59,*]

        map_2_7keV_tmp=make_array(60,30)
        map_2_7keV_tmp[0:15,*]=map_2_7keV[nflip:nflip+15,*]
        map_2_7keV_tmp[16:15+nflip,*]=map_2_7keV[0:nflip-1,*]
        map_2_7keV[59:0:-1,*]=map_2_7keV_tmp[0:59,*]

        map_4keV_tmp=make_array(60,30)
        map_4keV_tmp[0:15,*]=map_4keV[nflip:nflip+15,*]
        map_4keV_tmp[16:15+nflip,*]=map_4keV[0:nflip-1,*]
        map_4keV[59:0:-1,*]=map_4keV_tmp[0:59,*]
        
	; Add ribbon flux to map flux
	if band_num eq 0 then map1=map_1keV;/max(map_1keV)
	if band_num eq 1 then map1=map_2_7keV;/max(map_2_7keV)
	if band_num eq 2 then map1=map_4keV;/max(map_4keV)
endif

if gdf eq 1 then begin

        ; Read in GDF files
        map_0_7keV=make_array(30,60)
        map_1keV=make_array(30,60)
        map_1_7keV=make_array(30,60)
        map_2_7keV=make_array(30,60)
        map_4keV=make_array(30,60)

        file_0_7keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/f3-gdfMaps/hi09-gdf.txt'
        file_1keV='/home/kmarc/Downloads/RibbonSeparationPaperReleaseDataDirectories/f3-gdfMaps/hi3-gdf.txt'
        file_1_7keV='/home/kmarc/ENAMAPCODE_2016/RibbonSeparation/f3-gdfMaps/hi11-gdf.txt'
        file_2_7keV='/home/kmarc/Downloads/RibbonSeparationPaperReleaseDataDirectories/f3-gdfMaps/hi5-gdf.txt'
        file_4keV='/home/kmarc/Downloads/RibbonSeparationPaperReleaseDataDirectories/f3-gdfMaps/hi6_gdf.txt'

        map_0_7keV=read_table(file_0_7keV)
        map_1keV=read_table(file_1keV)
        map_1_7keV=read_table(file_1_7keV)
        map_2_7keV=read_table(file_2_7keV)
        map_4keV=read_table(file_4keV)

        ;map_0_7keV_tmp=map_0_7keV
        ;map_1keV_tmp=map_1keV
        ;map_1_7keV_tmp=map_1_7keV
        ;map_2_7keV_tmp=map_2_7keV
        ;map_4keV_tmp=map_4keV

        nflip=45.-1.

        map_0_7keV_tmp=make_array(60,30)
        map_0_7keV_tmp[0:15,*]=map_0_7keV[nflip:nflip+15,*]
        map_0_7keV_tmp[16:15+nflip,*]=map_0_7keV[0:nflip-1,*]
        map_0_7keV[59:0:-1,*]=map_0_7keV_tmp[0:59,*]

        map_1keV_tmp=make_array(60,30)
        map_1keV_tmp[0:15,*]=map_1keV[nflip:nflip+15,*]
        map_1keV_tmp[16:15+nflip,*]=map_1keV[0:nflip-1,*]
        map_1keV[59:0:-1,*]=map_1keV_tmp[0:59,*]

        map_1_7keV_tmp=make_array(60,30)
        map_1_7keV_tmp[0:15,*]=map_1_7keV[nflip:nflip+15,*]
        map_1_7keV_tmp[16:15+nflip,*]=map_1_7keV[0:nflip-1,*]
        map_1_7keV[59:0:-1,*]=map_1_7keV_tmp[0:59,*]

        map_2_7keV_tmp=make_array(60,30)
        map_2_7keV_tmp[0:15,*]=map_2_7keV[nflip:nflip+15,*]
        map_2_7keV_tmp[16:15+nflip,*]=map_2_7keV[0:nflip-1,*]
        map_2_7keV[59:0:-1,*]=map_2_7keV_tmp[0:59,*]

        map_4keV_tmp=make_array(60,30)
        map_4keV_tmp[0:15,*]=map_4keV[nflip:nflip+15,*]
        map_4keV_tmp[16:15+nflip,*]=map_4keV[0:nflip-1,*]
        map_4keV[59:0:-1,*]=map_4keV_tmp[0:59,*]

        ; Add ribbon flux to map flux
        if band_num eq 0 then map1=map_0_7keV
        if band_num eq 1 then map1=map_1keV
        if band_num eq 2 then map1=map_1_7keV
        if band_num eq 3 then map1=map_2_7keV
        if band_num eq 4 then map1=map_4keV

endif

if plot_limited eq 1 then begin
	plot_limits = [latmin,lonmin,latmax,lonmax]
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

if tail eq 1 then begin
	print, map_title
	print, 'nose', map1(58,14)  ;only for 6 degree pixels
	print, 'tail', map1(29,14)
	print, 'North pole', map1(58,29)
	print, 'South pole', map1(58,0)
	print, 'East Flank (right)', map1(8,14);map1(44,14)
	print, 'West Flank (left)', map1(38,14);map1(14,14)
	print, 'Voyager 1', map1(58,20)
	print, 'Voyager 2', map1(52,9)
endif else begin
	print, map_title
	print, 'nose', map1(29,14)  ;only for 6 degree pixels
	print, 'tail', map1(0,14)
	print, 'North pole', map1(29,29)
	print, 'South pole', map1(29,0)
	print, 'East Flank (right)', map1(39,14);map1(44,14)
	print, 'West Flank (left)', map1(19,14);map1(14,14)
	print, 'Voyager 1', map1(29,20)
	print, 'Voyager 2', map1(23,9)
endelse


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

window,0,retain=2
set_plot,'x'
device,true=24
device,decomposed=0
!p.background = 255.
!p.color = 0.
;!p.background = 255B
;!p.color = 0B
cgloadct,33,/silent,ncolors=254;loadct,3,/silent;
;Erase, Color=255;remove later
;include invisible title for spacing
;if plot_limited eq 1 then map_set,0,0,/HAMMER,/ISOTROPIC,/noerase, color=0, title = map_title,charsize = 2,limit = plot_limits else map_set,0,0,/HAMMER,/ISOTROPIC, color=0,charsize = 2 
if plot_limited eq 1 then map_set,0,0,/HAMMER,/ISOTROPIC,/noerase, color=0, title = map_title,charsize = 2,limit = plot_limits, /noborder else map_set,0,0,/HAMMER,/ISOTROPIC, color=0,charsize = 2
result = map_patch(map2,xstart=startx,ystart=starty,missing=255) ; missing=255 sets white background
;result(0,0) = 255.
thisposition=[0.0125,0.049,0.987,0.915]
;TV,result,startx,starty
cgimage,result,position=thisposition;,missing_color=cgcolor(annotateColor);,top=maxcolor;specifies the maximum color to plot to
;TVSCL,result,startx,starty,top=maxcolor ; top=maxcolor should be in there, but makes background black :(
;map_grid,latdel=map_dlat,londel=map_dlon,latlab = map_numlatlabels, lons = map_lons, lonnames =map_lonnames, color = 255, charsize =1.5,/label
;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Taken out by MZK on Dec 11 2015
;
;put on title on (can't have white lettering without box)
;!p.font = 1
xyouts,300,505,map_title,/device,color=0,charsize=3.0,charthick=2.
;!p.font  =-1
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

map_grid,color=0,glinethick=3,latdel=10,lats=[35];,45,55,65,75,85,95,105,115,125,135] ;Removed to get rid of grid lines

;plots voyager and nose positions, correct for nose rotated to 0 (should be 255)
noselong=-(255.-255.) ;0
taillong=-(75.-255.) ;180
voyager1long = -(253.-255.);-7.7
voyager2long=-(289.-255.);34
noselat = 5 ;just above grid lines
taillat= -5
voyager1lat= 35
voyager2lat=-32

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Taken out by MZK on Dec. 11 2015
;
if tail eq 1 then begin
	x= [noselong];
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
   	oplot, [x-.1,x,x+.1], [y-.1,y,y+.1],symsize=2, psym=8, color=0
    	xyouts,x+3.0,y+3.0, label, color=0,charsize=2
endif else begin
	for i=0, 2 do begin
	   oplot, [x[i]-.1,x[i],x[i]+.1], [y[i]-.1,y[i],y[i]+.1],symsize=2, psym=8, color=0
	    xyouts,x[i]+3.0,y[i]+3.0, label[i], color=0,charsize=2
	endfor
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

!X.TICKFORMAT = ''

;make colorbar
if plot_log ne 1 then begin
	if units eq 'Counts/s' then cgcolorbar, /vertical, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(f10.3)', position =[0.69,0.93,0.99,0.99], charthick=2.,charsize =2 else cgcolorbar, /vertical, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(f10.1)', position =[0.79,0.935,0.99,0.995],charthick=2, charsize =2 ;( colorbar.pro at dfanning.com)
endif else begin
    ;cb_ticknames =string((maxset/minset)^(findgen(5)/4.)*minset,format='(I3)')
    cb_ticknames =string((maxset-minset)*(findgen(5)/4.)^2.+minset,format='(I3)')

	if units eq 'Counts/s' then colorbar, /vertical, bottom = 0, divisions = 4., minor =0., ticknames= cb_ticknames, format = '(f10.3)', position =[0.69,0.93,0.99,0.99], charsize =1.5 else colorbar, /vertical, bottom = 0, divisions = 4., minor =0., ticknames= cb_ticknames, format = '(f10.1)', position =[0.69,0.935,0.99,0.995], charsize =1.5  ;( colorbar.pro at dfanning.com)

endelse


;save
image = TVRD(0,0,!D.x_size,!D.y_size,true=3)
WRITE_JPEG, filename+'.jpg', image, QUALITY=100, TRUE=3

end
;-------------------------------------------
FUNCTION ALOGNUM, x
common log_def, log_num
      RETURN, ALOG(x) / ALOG(log_num)
END

