pro get_limits

common rlimits, r_limits_lower, r_limits_upper, boundaries_arr
common magnetic_field, ophermagnetic
common grid_params, dtheta, dphi, dr, nphi, ntheta, phi_first, phi_last,  theta_first,  theta_last, ri_opher, rf_opher,  r_lower,r_upper,r_inner,moscow_regs
common secondary_path_names,secondary_path, opher_withneutrals_secondary_file, opher_withneutrals_secondarynh_file, opher_withneutrals_secondarystream_file,input_path,output_path
common grid_opher,gridopher
common radial_limits, radial_range

nr=gridopher.r.num; number of r bins

r_limits_lower = fltarr(ntheta,nphi) & r_limits_upper = fltarr(ntheta,nphi)

boundaries_arr = intarr(ntheta,nphi,nr)
boundaries_arr0 = intarr(ntheta,nphi,nr)
other=dblarr(ntheta,nphi,nr)
ophermagnetic = {magnetic, Bx:other, By:other, Bz:other}

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; If using boundaries.txt as derived from Tepclot, will have indices outputted as floats
; instead of integers. This accomodates for that fact.
; - MZK 2/10/16
;
;i_temp = 0 & j_temp = 0 &  k_temp=0 & theta_temp=0. & phi_temp = 0. & r_temp = 0. & region_temp = 0. & Bx_temp = 0. & By_temp = 0. & Bz_temp = 0.

i_temp = 0. & j_temp = 0. &  k_temp=0. & theta_temp=0. & phi_temp = 0. & r_temp = 0. & region_temp = 0. & Bx_temp = 0. & By_temp = 0. & Bz_temp = 0.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

close,1
openr,1,secondary_path+input_path+'boundaries.txt'
;for i=0, ntheta-1 do begin
;    theta=theta_first+i*dtheta
;    for j=0, nphi-1 do begin
;        phi=phi_first+j*dphi
;        for k=0, nr-1 do begin
;            r=ri_opher+k*dr

for i=0, ntheta-1 do begin ; For Tecplot outputs - MZK
     for k=0, nr-1 do begin
	 for j=0, nphi-1 do begin
	    phi=phi_first+j*dphi 
	    theta=theta_first+i*dtheta 
	    r=ri_opher+k*dr  
            readf,1,region_temp,i_temp,j_temp,k_temp,theta_temp,phi_temp,r_temp,Bx_temp, By_temp, Bz_temp		
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; In For loop, need to accomodate for indices being floats from Tecplot 
; -MZK 2/10/16
;
;            if i_temp ne i or j_temp ne j or k_temp ne k or theta_temp ne theta or phi_temp ne phi or r_temp ne r then stop 
;            r_temp=r_temp+30 ; Needed this for single ion. ******LOOK INTO WHY*****

	    if theta_temp-floor(theta_temp) lt 0.5 then t=floor(theta_temp)+0.
	    if theta_temp-floor(theta_temp) ge 0.5 then t=ceil(theta_temp)+0.
            if phi_temp-floor(phi_temp) lt 0.5 then p=floor(phi_temp)+0.
            if phi_temp-floor(phi_temp) ge 0.5 then p=ceil(phi_temp)+0.
	    if r_temp-floor(r_temp) lt 0.5 then rad=floor(r_temp)+0.
	    if r_temp-floor(r_temp) ge 0.5 then rad=ceil(r_temp)+0.
;	    t=floor(theta_temp)+0. ; added to make boundaries file directly from interpolated 3D work - MZK
;	    p=floor(phi_temp)+0.
;	    t=ceil(theta_temp)+0. ; Rounds phi up since slighly off from actual value (from below)
;           p=ceil(phi_temp)+0. ; Rounds phi up since slighly off from actual value (from below)
	    a=i+0. ; turns i into float for if statement
	    b=j+0. ; turns j into float for if statement
	    c=k+0. ; turns k into float for if statement
;	    if i_temp ne i or j_temp ne j or k_temp ne k or theta_temp ne theta or phi_temp ne phi or r_temp ne r then stop
;	    if i_temp ne a or j_temp ne b or k_temp ne c or t ne theta or p ne phi or rad ne r then stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 	    boundaries_arr0[i,j,k] = region_temp
   	    if moscow_regs eq 1 then begin
		if boundaries_arr0[i,j,k] eq 3 then boundaries_arr[i,j,k]=2
		if boundaries_arr0[i,j,k] eq 4 then boundaries_arr[i,j,k]=2
		if boundaries_arr0[i,j,k] eq 1 then boundaries_arr[i,j,k]=4
		if boundaries_arr0[i,j,k] eq 2 then boundaries_arr[i,j,k]=3
	    endif else boundaries_arr[i,j,k] = boundaries_arr0[i,j,k]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; USE FOR ALL FILES THAT ARE FROM TECPLOT OUTPUTS
;
;if boundaries_arr[i,j,k] ne 1. or boundaries_arr[i,j,k] ne 2. or boundaries_arr[i,j,k] ne 3. or boundaries_arr[i,j,k] ne 4. or boundaries_arr[i,j,k] ne 5. then boundaries_arr[i,j,k]=5.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


            ophermagnetic.Bx[i,j,k]= Bx_temp &  ophermagnetic.By[i,j,k]= By_temp & ophermagnetic.Bz[i,j,k]= Bz_temp ; in Gauss
        endfor
    endfor
endfor
close,1

for i=0, ntheta-1 do begin
    for j=0, nphi-1 do begin
		flag_found_begin = 0.
		flag_found_end = 0.
		for k=0, nr-1 do begin
		
		;TS to HP
		if radial_range eq 0 then begin			
	
                        if flag_found_begin eq 0. then if boundaries_arr[i,j,k] eq 3 then begin
				if ri_opher+k*dr gt 60. then begin  ;check for bad grid
					r_limits_lower(i,j) =  k;ri_opher+k*dr ; EXPLORE DEFINITION OF THIS FOR makecountratefiles-MZK
					flag_found_begin  = 1.	
				endif		
			endif
                       if flag_found_end eq 0. then if boundaries_arr[i,j,k] eq 2 then begin
				r_limits_upper(i,j) =  k-1;ri_opher+(k-1)*dr
				flag_found_end = 1.		
                        endif
                       if flag_found_end eq 1. and boundaries_arr[i,j,k] eq 3 then begin
                           flag_found_end = 0.
;                           print, 'changing flag', i,j,k
                          
                       endif
		endif


                ;HP to BS (800)
		if radial_range eq 1 then begin		
			if flag_found_begin eq 0. then if boundaries_arr[i,j,k] eq 2 then begin
				if ri_opher+k*dr gt 80. then begin  ;check for bad grid
					r_limits_lower(i,j) =  k;ri_opher+k*dr
					flag_found_begin  = 1.	
                    		endif		
			endif    
                     	if flag_found_begin eq 1. and boundaries_arr[i,j,k] eq 3 then begin
                          flag_found_begin = 0.
;                          print, 'changing flag', i,j,k
                          
                       endif
			if flag_found_end eq 0. then if boundaries_arr[i,j,k] eq 1 then begin
				r_limits_upper(i,j) =  k-1;ri_opher+(k-1)*dr
				flag_found_end = 1.		
			endif	 
	
		endif

		if radial_range eq 2 then begin
		;TS to 800
        		   if flag_found_begin eq 0. then if boundaries_arr[i,j,k] eq 3 then begin
				if ri_opher+k*dr gt 60. then begin  ;check for bad grid
					r_limits_lower(i,j) =  k;ri_opher+k*dr
					flag_found_begin  = 1.	
				endif		
			   endif
		endif

		endfor
		if flag_found_begin eq 0 then r_limits_lower(i,j) = nr-1		
		if flag_found_end eq 0 then r_limits_upper(i,j) = nr-1
               
;                print, ri_opher+r_limits_lower(i,j)*dr, ri_opher+r_limits_upper(i,j)*dr
      	endfor

endfor

return
end
