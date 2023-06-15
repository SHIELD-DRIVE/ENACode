pro plotter_eps, map,filename, map_title, units, maxset, minset, band_num
common plot_params,latmin,latmax,lonmin,lonmax,map_dlat,map_dlon,map_numlatlabels,map_lons,map_lonnames, onerange_flux_min, onerange_flux_max, onerange_cr_min, onerange_cr_max
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher
common plot_flags, onerange_flag,manual_onerange_flag, plot_log,plot_limited, manual_lim, manual_maxset, manual_minset, eps_plot, plot_scl,tail,legend_style
common model_names,models,num_models, colde

; flips array upside down (so that it plots 0 theta at the top of the map)
map1=fltarr(nphi,ntheta)  ;flipped for mappings nomenclature

print, "ENERGY BAND: ", map_title

print, "Nose flux", map[14,30]
print, "Tail flux", map[14,0]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 

map1tmp=map*plot_scl

for i=0,ntheta-1 do begin
    for j=0, nphi-1 do begin
        map1[j,(ntheta-1-i)]=map1tmp[i,j]
    end
end

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
   for j=0, nflip-1 do begin
      map1[j,*]=map4atmp[nflip+j,*]
      map1[nflip+j,*]=map4atmp[j,*]
   endfor

map4btmp=map1                    ; changing outside-in to inside-out for eps
for j=nphi-1, 0, -1 do map1[nphi-1-j,*]=map4btmp[j,*]

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if plot_limited eq 1 then begin
        plot_limits =[latmin,lonmin,latmax,lonmax]
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
    ;sqrt scheme
    if manual_lim eq 1 then maxcolor=max(map2) ; if using manual plot limits - MZK
    if manual_lim ne 1 then maxcolor = 253 ; if using default plot limits - MZK
endif

;window,0,retain=2 ; gets plot formatting correct without jpeg

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
if tail eq 1 then begin
	if plot_limited eq 1 then cgmap_set,-5,79,/HAMMER,/ISOTROPIC,/noerase, color=0, charsize = 2.5,position=[0.0148, 0.214, 0.983, 0.942],limit = plot_limits, /noborder else map_set,-5,79,/HAMMER,/ISOTROPIC, color=0,charsize = 2.5, position=[0.0148, 0.214, 0.983, 0.942],/noborder, reverse=1
endif else begin
	if plot_limited eq 1 then cgmap_set,5,-101,/HAMMER,/ISOTROPIC,/noerase, color=0, charsize = 2.5,position=[0.0148, 0.214, 0.983, 0.942],limit = plot_limits, /noborder else map_set,5,-101,/HAMMER,/ISOTROPIC, color=0,charsize = 2.5, position=[0.0148, 0.214, 0.983, 0.942],/noborder, reverse=1
endelse
result = map_image(map2,startx,starty,missing=255,compress=1);,/bilinear) ; include bilinear for cell interpolation
thisposition=[0.02, 0.212, 0.97, 0.94] 
cgimage,result,position=thisposition

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

;make colorbar
if plot_log ne 1 then begin
	if legend_style eq 0 then if units eq 'Counts/s' then cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(f10.3)', position =[0.69,0.93,0.99,0.99], charthick=4, charsize =3 else cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(I4)', position =[0.1,0.17,0.9,0.2], charthick=8, charsize =3. ;( colorbar.pro at dfanning.com)
       if legend_style eq 1 then if units eq 'Counts/s' then cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(f10.3)', position =[0.69,0.93,0.99,0.99], charthick=4, charsize =3 else cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(f10.2)', position =[0.1,0.17,0.9,0.2], charthick=8, charsize =3. ;( colorbar.pro at dfanning.com)
       if legend_style eq 2 then if units eq 'Counts/s' then cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(f10.3)', position =[0.69,0.93,0.99,0.99], charthick=4, charsize =3 else cgcolorbar, bottom = 0, divisions = 4., minor =0., range = [minset,maxset], format = '(e9.1)', position =[0.1,0.17,0.9,0.2], charthick=8, charsize =3. ;( colorbar.pro at dfanning.com)
endif else begin
    cb_ticknames =string((maxset-minset)*(findgen(5)/4.)^2.+minset,format='(I3)')

	if units eq 'Counts/s' then cgcolorbar, bottom = 0, divisions = 4., minor =0., ticknames= cb_ticknames, format = '(f10.3)', position =[0.69,0.93,0.99,0.99], charthick=2, charsize =2 else cgcolorbar, bottom = 0, divisions = 4., minor =0., ticknames= cb_ticknames, format = '(I4)', position =[0.75,0.945,0.99,0.995], charthick=3, charsize =2.5  ;( colorbar.pro at dfanning.com)

endelse

ps_off

print, ""

end
;-------------------------------------------
FUNCTION ALOGNUM, x
common log_def, log_num
      RETURN, ALOG(x) / ALOG(log_num)
END

