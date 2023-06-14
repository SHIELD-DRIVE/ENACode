pro read_stream, lc, cool, vx, vy, vz, a, np_stream, tp_stream, vr_stream, vt_stream
  
  common secondary_path_names,secondary_path, opher_secondary_file, heerikhuisen_secondary_file, heerikhuisen_secondarynh_file, opher_withneutrals_secondary_file, opher_withneutrals_secondarynh_file, opher_twoplasmas_secondaryp1_file,opher_twoplasmas_secondaryp2_file, opher_twoplasmas_secondarynh_file, opher_withneutrals_secondarystream_file,opher_withneutrals_secondary_file1, opher_withneutrals_secondarynh_file1,opher_withneutrals_secondary_file2, opher_withneutrals_secondarynh_file2,opher_withneutrals_secondaryinterp_file,input_path,output_path
  common constants_cgs, mp,kb,kev_erg, AU_cm
  common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen,r_lower,r_upper,r_inner,moscow_regs
  common grid_opher,gridopher
  common rlimits, r_limits_lower, r_limits_upper, boundaries_arr
  common model_names, models, num_models, pui_model, vasyliunas, ribbon,gdf, hinterp, ntr_frac, nref_frac, ttr_frac, tref_frac,colde,esw_frac,etr_frac,eref_frac,moscow,moscow_temp
  
  nr=gridopher.r.num; number of r bins
  sv=7.7437451e-08 ; vena(1 keV)*sigma(1 keV)
  v1kev=4.608e7 
  
  file=secondary_path+input_path+opher_withneutrals_secondarystream_file

  data=read_ascii(file,data_start=25)

  x=data.field01(0,*)            ; au
  y=data.field01(1,*) ; au
  z=data.field01(2,*) ; au
  ux=data.field01(3,*)*1.e5 ; cm/s
  uy=data.field01(4,*)*1.e5 ; cm/s
  uz=data.field01(5,*)*1.e5 ; cm/s
  nh=data.field01(6,*)     ; cc
  un_x=data.field01(7,*)*1.e5 ; cm/s
  un_y=data.field01(8,*)*1.e5 ; cm/s
  un_z=data.field01(9,*)*1.e5 ; cm/s
  np=data.field01(10,*) ; cc
  if moscow_temp eq 1 then tp=data.field01(11,*) else tp=data.field01(11,*)/2. ; K
  vr=data.field01(12,*)*1.e5 ; cm/s
  vtheta=data.field01(13,*)*1.e5 ; cm/s
  vphi=data.field01(14,*)*1.e5 ; cm/s
  vt=data.field01(15,*)*1.e5 ; cm/s
  theta=data.field01(16,*)      ; deg
  phi=data.field01(17,*) ; deg
  r=data.field01(18,*)   ; au
  if moscow_regs eq 1 then reg=data.field01(19,*)

  nan=where(np ne np)
  num=n_elements(np)
  ext=dblarr(num-1)
  nh_sum=dblarr(num-1)
  dx_sum=dblarr(num-1)
  a0=dblarr(num-1)
  np0_stream=dblarr(num-1)
  tp0_stream=dblarr(num-1)
  vr0_stream=dblarr(num-1)
  vt0_stream=dblarr(num-1)
  ptr_ct=dblarr(ntheta,nphi,nr)
  ptra_ct=dblarr(ntheta,nphi,nr)
  ptrn_ct=dblarr(ntheta,nphi,nr)
  btest=dblarr(ntheta,nphi,nr)

;fixing phi
  phi=atan(y,x)*57.2958
  for q=0, num-1 do if phi[q] lt 0 then phi[q]=360.-abs(phi[q])
  
  v=sqrt(vr^2.+vt^2.)
  cs=sqrt(5*(2*tp)*kb/(3*mp)) ; multiply by 2 to account for electrons
  M2=(v^2.)/(cs^2.)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; If read in energy, then use this
; If don't read in energy, comment out
;
;  
  ; pop 1 radial and tangential velocities
  vn_r=(un_x*x+un_y*y+un_z*z)/r
  vn_t=sqrt((un_x^2+un_y^2+un_z^2)-vn_r^2)

  ; Making arrays for different neutrals pops
  w=dblarr(num);,4)
  vrel=dblarr(num);,4)
  erel=dblarr(num);,4)
  sig1=dblarr(num);,4)
  sig2=dblarr(num);,4)
  sig=dblarr(num);,4)
  
  vth=sqrt(kb*tp/mp)

  ; Using parent proton velocity
  w(*)=SQRT((un_x-ux)^2.+(un_y-uy)^2.+(un_z-uz)^2.)/vth

  ; Relative velocity between neutrals and ions
     vrel=vth*(exp(-w^2.)/sqrt(!pi)+(w+1./(2.*w))*erf(w))

  ; Cross section
     a1=4.15
     a2=0.531
     a3=67.3

     erel=0.5*mp*vrel^2.
     erel=erel/kev_erg
  
     sig1=(a1-a2*alog(erel))^2.
     sig2=(1-exp(-a3/erel))^(4.5)
     sig=sig1
     sig=sig*10D^(-16D)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dr0=r[1:num-1]-r[0:num-2]
  dr0=dr0*au_cm ; convert from au to cm
  dtheta0=theta[1:num-1]-theta[0:num-2]
  dtheta0=dtheta0*(!pi/180.) ; convert from deg to rad
  dphi0=phi[1:num-1]-phi[0:num-2]
  dphi0=dtheta0*(!pi/180.) ; convert from deg to rad
  ds=sqrt(dr0^2.+r[1:num-1]^2.*dtheta0^2.+r[1:num-1]^2.*sin(theta[1:num-1]*(!pi/180.))^2.*dphi0^2.) ; spherical polar line element

 tau=-(nh[0:num-2]/v[1:num-1])*ds[0:num-2]

 nu1au=8.e-8                    ; 1/s
 b=0.08
 
if moscow_regs eq 1 then begin
  for i=1,num-2 do begin
     if reg[i-1] eq 2 and reg[i] lt 2 and reg[i] gt 1 then reg[i]=2.
     if reg[i] ge 1 and reg[i] lt 2 and r[i] ge r_lower and r[i] lt 200. then begin
         dx_sum[i]=(r[i+1]-r[i])+dx_sum[i-1] ; determining factors for weighted avg. of neutrals
         nh_sum[i]=nh[i]*(r[i+1]-r[i])+nh_sum[i-1]
     endif
     ;if reg[i] ge 1 and reg[i] lt 2 and r[i] lt r_lower+3.*dr and M2[i-1] eq M2[i-1] or i eq 1 then begin
     if reg[i] ge 1 and reg[i] lt 2 and r[i] lt r_lower+2.*dr and M2[i-1] eq M2[i-1] or i eq 1 then begin
        u1au=v[i]
        n1au=np[i]*r[i]*r[i]
        sig_cx_sw=sig[i]
     endif
     if reg[i] ge 1 and reg[i] lt 2 and reg[i+1] eq 2  and r[i] lt 200. then begin
        nh_avg=nh_sum[i]/dx_sum[i] ; weighted average of neutrals in supersonic SW
        r_ts=r[i]*au_cm
        a0[i]=r_ts*nh_avg*(nu1au+sig_cx_sw*u1au*n1au)/(u1au*n1au)
        if a0[i] gt 1 then stop ; 316611
        if a0[i] ne a0[i] then stop
        ;if a0[i] eq a0[i] then print,i, a0[i], r_ts/au_cm, theta[i], phi[i], u1au, n1au, nh_avg
     endif else if reg[i] eq 2 and r[i] eq r[i] and r[i-1] eq r[i-1] then a0[i]=a0[i-1]
  endfor

  nan=0./0.
; Ensuring only considering points starting from TS
  for i=1, num-2 do begin
     if reg[i] ge 1 and reg[i] lt 2 and r[i] lt 200. then tau[i]=0
     if (r[i] eq r[i]) and (r[i-1] eq r[i-1]) then ext[i]=tau[i]+ext[i-1] $
     else ext[i]=0
  endfor

endif else begin
  for i=1,num-2 do begin
     if M2[i] gt 1.0 and r[i] ge r_lower and r[i] lt 130. then begin
         dx_sum[i]=(r[i+1]-r[i])+dx_sum[i-1] ; determining factors for weighted avg. of neutrals
         nh_sum[i]=nh[i]*(r[i+1]-r[i])+nh_sum[i-1]
     endif
     if M2[i] gt 1.0 and r[i] lt r_lower+3.*dr and M2[i-1] eq M2[i-1] or i eq 1 then begin
        u1au=v[i]
        n1au=np[i]*r[i]*r[i]
        sig_cx_sw=sig[i]
     endif
     if M2[i] gt 1.0 and M2[i+1] lt 1.0 and r[i] lt 130 then begin
        nh_avg=nh_sum[i]/dx_sum[i] ; weighted average of neutrals in supersonic SW
        r_ts=r[i]*au_cm
        if M2[i] gt 1.0 then a0[i]=r_ts*nh_avg*(nu1au+sig_cx_sw*u1au*n1au)/(u1au*n1au)
        if a0[i] gt 1 then stop
     endif else if M2[i] lt 1.0 and r[i] eq r[i] and r[i-1] eq r[i-1] then a0[i]=a0[i-1]
  endfor  

  nan=0./0.
; Ensuring only considering points starting from TS
  for i=1, num-2 do begin
     if M2[i] gt 0.81 then tau[i]=0
     if (r[i] eq r[i]) and (r[i-1] eq r[i-1]) then ext[i]=tau[i]+ext[i-1] $
     else ext[i]=0
  endfor
endelse

  ;for i=1, num-2 do begin ; For following density along streamline for extinction
  ;   if (M2[i] lt 0.81) and (M2[i-1] ge 0.81) and (r[i] lt 200.) then begin
  ;      zm0=where(np[i:i+15]*exp(ext[i:i+15]*sv) eq max(np[i:i+15]*exp(ext[i:i+15]*sv)))
  ;      zm=i+zm0[0]
  ;      i=zm
  ;   endif
  ;   if M2[i] lt 0.81 then np0_stream[i]=np[zm]
  ;   if M2[i] lt 0.81 then tp0_stream[i]=tp[zm]
  ;   if M2[i] lt 0.81 then vr0_stream[i]=vr[zm]
  ;   if M2[i] lt 0.81 then vt0_stream[i]=vt[zm]
  ;endfor
  
  if num lt 2000 then streamline_test, ext, r, nh, v, vr, vt, np, tp
;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; For ApJ paper
;enum=13
;  energy=[0.1,0.7,1.1,1.7,2.7,4.29,10,20,30,40,50,75,100]
;  
;  vena=energy_velocity(energy)
;
;  a1=4.15
;  a2=0.531
;  a3=67.3
;
;  sig1=(a1-a2*alog(energy))^2
;  sig2=(1-exp(-a3/energy))^(4.5)
;  sigma=sig1*sig2
;  sigma=sigma*10D^(-16D)
;
;  ext_arr=make_array(enum,num-1)
;  for i=0,enum-1 do ext_arr(i,*)=ext[*]*vena[i]*sigma[i]
;
; ; 6 deg: 50642-->51385 streamline w/ TS @ theta=71.84 deg and phi=6.01880 deg
; ; 3 deg: 190094-->190843 streamline w/ TS @ theta=68.24 deg and phi=5.99 
;
;  lc=make_array(enum)
;  for j=0,enum-1 do begin 
;	depth=where(abs(1/exp(1)+ext_arr[j,190094:190843]) eq min(abs(1/exp(1)+ext_arr[j,190094:190843])),dct)
;	lc[j]=min(depth)
;  endfor
;
;  lc=r(190094+lc)
;
;  set_plot, 'ps'
;  device, filename="cooling_length_apj.ps", /color, /encapsulated
;  device, /helvetica, /bold
;  !p.font=0
;  loadct, 39
;  plot, [0.1,100], [100,1000],/NODATA, ytitle='Cooling Length [AU]',xtitle='Energy [keV]', xstyle=1,ystyle=1, xthick=3,ythick=3,ylog=1,xlog=1,Position=[0.16, 0.15, 0.84, 0.85],charsize=1.15
;  oplot, energy, lc, thick=5, psym=-4
;  device,/close
;
;  set_plot, 'ps'
;  device, filename="streamline_apj.ps", /color, /encapsulated
;  device, /helvetica, /bold
;  !p.font=0
;  loadct, 39
;  plot, [0,1500], [0,400],/NODATA, ytitle='Speed [km/s]',xtitle='Radius Along Streamline [AU]', xstyle=1,ystyle=1, xthick=3,ythick=3,ylog=0,xlog=0,Position=[0.16, 0.15, 0.84, 0.85],charsize=1.15
;  oplot, r[190094:190843], v[190094:190843]/1.e5, thick=5
;  device,/close 
;
;  stop
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dist=where(ext gt -1.01 and ext lt -0.99,ct)
  r_in=round((r-r_lower)/dr) ; indices for radius
  t_in=round((theta-2.9)/dtheta) ; indices for theta
  p_in=round((phi-2.9)/dphi)     ; indices for phi
;  t_in=round(theta/dtheta)
;  p_in=round(phi/dphi)

 rlower=(r_lower-r_inner)/dr; 23
 rupper=(r_upper-r_inner)/dr

for i=0,10 do begin
  for t=0,ntheta-1 do begin
     for p=0,nphi-1 do begin
        if i eq 0 then begin ; putting extinction values into 3D array of theta, phi and r
           print, string((nphi*t+p)/(nphi*ntheta)*100.)+'%'
           for rad=rlower, rupper do begin
             if boundaries_arr[t,p,rad] eq 3. or boundaries_arr[t,p,rad] eq 7 then begin
                 ptr0=where(p_in eq p and t_in eq t and r_in eq rad,ptr0ct) ; Actually 'ptr-1' based on the way tau and ext is calculated
		 if ptr0ct gt 0 then ptr_test=where(ext[ptr0] eq ext[ptr0],ptrct) else ptrct=0
                 if ptrct gt 0 then ptr=ptr0[ptr_test] else ptr=-1
		 if ptrct gt 0 then ptr_ct[t,p,rad]=ptrct
                 if ptrct gt 0 then cool[t,p,rad]=mean(ext[ptr])
                 if cool[t,p,rad] ne cool[t,p,rad] and t eq 0 and p eq 0 then cool[t,p,rad]=0.
                 if ptr(0) ge 0 then vx[t,p,rad]=mean(ux[ptr]) else vx[t,p,rad]=0. ; for printing out streamlines
                 if ptr(0) ge 0 then vy[t,p,rad]=mean(uy[ptr]) else vy[t,p,rad]=0.
                 if ptr(0) ge 0 then vz[t,p,rad]=mean(uz[ptr]) else vz[t,p,rad]=0.

		 if ptr0ct gt 0 then ptra_test=where(a0[ptr0] eq a0[ptr0] and a0[ptr0] ne 0,ptract) else ptract=0
                 if ptract gt 0 then ptra=ptr0[ptra_test] else ptra=-1
		 if ptract gt 0 then ptra_ct[t,p,rad]=ptract
                 if ptra(0) ge 0 then a[t,p,rad]=mean(a0[ptra])

                 if ptr0ct gt 0 then ptrn_test=where(np0_stream[ptr0] eq np0_stream[ptr0] and np0_stream[ptr0] ne 0,ptrnct) else ptrnct=0
                 if ptrnct gt 0 then ptrn=ptr0[ptrn_test] else ptrn=-1
		 if ptrnct gt 0 then ptrn_ct[t,p,rad]=ptrnct
                 if ptrn(0) ge 0 then np_stream[t,p,rad]=mean(np0_stream[ptrn])
                 if ptrn(0) ge 0 then tp_stream[t,p,rad]=mean(tp0_stream[ptrn])
                 if ptrn(0) ge 0 then vr_stream[t,p,rad]=mean(vr0_stream[ptrn])
                 if ptrn(0) ge 0 then vt_stream[t,p,rad]=mean(vt0_stream[ptrn])
              endif else cool[t,p,rad]=0.                 
           endfor
        endif

        if i eq 1 then begin    ; first wave of interpolating zero extinction values
           for rad=rlower, rupper do begin
              if boundaries_arr[t,p,rad] eq 3. then begin
                 n=0.
                 na=0.
                 nd=0.
                 
                 if cool[t,p,rad] eq 0 and t ne 0 and p ne 0 and t ne ntheta-1 and p ne nphi-1 then begin ; interpolating for 3 deg < theta < 357 deg and 3 deg < phi < 357 deg
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                          if cool[t+q,p+u,rad] ne 0 and ptr_ct[t+q,p+u,rad] gt 0 and abs(q)+abs(u) ne 0 and cool[t+q,p+u,rad] eq cool[t+q,p+u,rad] then begin
                             cool[t,p,rad]=cool[t+q,p+u,rad]+cool[t,p,rad] ; using all surrounding theta, phi values for interpolation
                             n=n+1.                                        ; if nonzero extinction, adds to index for averaging
                          endif
                       endfor
                    endfor
                    if n ne 0 then begin
                       cool[t,p,rad]=cool[t,p,rad]/n
                    endif
                 endif

                 if a[t,p,rad] eq 0 and t ne 0 and p ne 0 and t ne ntheta-1 and p ne nphi-1 then begin ; interpolating for 3 deg < theta < 357 deg and 3 deg < phi < 357 deg
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if a[t+q,p+u,rad+g] ne 0 and ptra_ct[t+q,p+u,rad+g] gt 0 and abs(q)+abs(u)+abs(g) ne 0 and a[t+q,p+u,rad+g] eq a[t+q,p+u,rad+g] then begin
                                a[t,p,rad]=a[t+q,p+u,rad+g]+a[t,p,rad]
                                na=na+1. ; if nonzero extinction, adds to index for averaging
                             endif
                          endfor
                       endfor
                    endfor
                    if na ne 0 then begin
		       a[t,p,rad]=a[t,p,rad]/na
                    endif
                 endif

                 if np_stream[t,p,rad] eq 0 and t ne 0 and p ne 0 and t ne ntheta-1 and p ne nphi-1 then begin ; interpolating for 3 deg < theta < 357 deg and 3 deg < phi < 357 deg
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if np_stream[t+q,p+u,rad+g] ne 0 and ptrn_ct[t+q,p+u,rad+g] gt 0 and abs(q)+abs(u)+abs(g) ne 0 and np_stream[t+q,p+u,rad+g] eq np_stream[t+q,p+u,rad+g] then begin
                                np_stream[t,p,rad]=np_stream[t+q,p+u,rad+g]+np_stream[t,p,rad]
                                tp_stream[t,p,rad]=tp_stream[t+q,p+u,rad+g]+tp_stream[t,p,rad]
                                vr_stream[t,p,rad]=vr_stream[t+q,p+u,rad+g]+vr_stream[t,p,rad]
                                vt_stream[t,p,rad]=vt_stream[t+q,p+u,rad+g]+vt_stream[t,p,rad]
                                nd=nd+1. ; if nonzero extinction, adds to index for averaging
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
		       np_stream[t,p,rad]=np_stream[t,p,rad]/nd
		       tp_stream[t,p,rad]=tp_stream[t,p,rad]/nd
		       vr_stream[t,p,rad]=vr_stream[t,p,rad]/nd
		       vt_stream[t,p,rad]=vt_stream[t,p,rad]/nd
                    endif
                 endif
              endif else cool[t,p,rad]=0. 
           endfor
        endif  
        
        if i ge 2 then begin
           for rad=rlower, rupper do begin
              if boundaries_arr[t,p,rad] eq 3. then begin
                 n=0.
                 na=0.
                 nd=0.
                 qp=0.
                               
                 if cool[t,p,rad] eq 0 and t ne 0 and p eq 0 and t ne ntheta-1 then begin ; interpolating for phi = 3 deg and 3 deg < theta < 357 deg (extinction)
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                          if u lt 0 and cool[t+q,nphi-1,rad] ne 0 and cool[t+q,nphi-1,rad] eq cool[t+q,nphi-1,rad] then begin
                             cool[t,p,rad]=cool[t+q,nphi-1,rad]+cool[t,p,rad]
                             n=n+1
                          endif else if u ge 0 then begin
                             if cool[t+q,p+u,rad] ne 0 and abs(q)+abs(u) ne 0 and cool[t+q,p+u,rad] eq cool[t+q,p+u,rad] then begin
                                cool[t,p,rad]=cool[t+q,p+u,rad]+cool[t,p,rad]
                                n=n+1.
                             endif
                          endif
                       endfor
                    endfor
                    if n ne 0 then begin
                       cool[t,p,rad]=cool[t,p,rad]/n
                    endif
                 endif

                 if a[t,p,rad] eq 0 and t ne 0 and p eq 0 and t ne ntheta-1 then begin ; interpolating for phi = 3 deg and 3 deg < theta < 357 deg (Density)
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if u lt 0 and a[t+q,nphi-1,rad+g] ne 0 and a[t+q,nphi-1,rad+g] eq a[t+q,nphi-1,rad+g] then begin
                                a[t,p,rad]=a[t+q,nphi-1,rad+g]+a[t,p,rad]
                                na=na+1
                             endif else if u ge 0 then begin
                                if a[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and a[t+q,p+u,rad+g] eq a[t+q,p+u,rad+g] then begin
                                   a[t,p,rad]=a[t+q,p+u,rad+g]+a[t,p,rad]
                                   na=na+1.
                                endif
                             endif
                          endfor
                       endfor
                    endfor
                    if na ne 0 then begin
                       a[t,p,rad]=a[t,p,rad]/na
                    endif
                 endif

                 if np_stream[t,p,rad] eq 0 and t ne 0 and p eq 0 and t ne ntheta-1 then begin ; interpolating for phi = 3 deg and 3 deg < theta < 357 deg (Density)
                    for q=-1,1 do begin
                       for u=0,1 do begin
                          for g=-1,1 do begin
                             if np_stream[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and np_stream[t+q,p+u,rad+g] eq np_stream[t+q,p+u,rad+g] then begin
                                np_stream[t,p,rad]=np_stream[t+q,p+u,rad+g]+np_stream[t,p,rad]
                                tp_stream[t,p,rad]=tp_stream[t+q,p+u,rad+g]+tp_stream[t,p,rad]
                                vr_stream[t,p,rad]=vr_stream[t+q,p+u,rad+g]+vr_stream[t,p,rad]
                                vt_stream[t,p,rad]=vt_stream[t+q,p+u,rad+g]+vt_stream[t,p,rad]
                                nd=nd+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       np_stream[t,p,rad]=np_stream[t,p,rad]/nd
                       tp_stream[t,p,rad]=tp_stream[t,p,rad]/nd
                       vr_stream[t,p,rad]=vr_stream[t,p,rad]/nd
                       vt_stream[t,p,rad]=vt_stream[t,p,rad]/nd
                    endif
                 endif
                 
                 if cool[t,p,rad] eq 0 and t eq ntheta-1 and p eq nphi-1 then begin ; interpolating for theta = 357 deg and phi = 357 deg (Extinction)
                    for q=-1,0 do begin
                       for u=-1,1 do begin
                          if u gt 0 and cool[t+q,0,rad] ne 0 and cool[t+q,0,rad] eq cool[t+q,0,rad] then begin
                             cool[t,p,rad]=cool[t+q,0,rad]+cool[t,p,rad]
                             n=n+1
                          endif else if u le 0 then begin
                             if cool[t+q,p+u,rad] ne 0 and abs(q)+abs(u) ne 0 and cool[t+q,p+u,rad] eq cool[t+q,p+u,rad] then begin
                                cool[t,p,rad]=cool[t+q,p+u,rad]+cool[t,p,rad]
                                n=n+1.
                             endif
                          endif
                       endfor
                    endfor
                    if n ne 0 then begin
                       cool[t,p,rad]=cool[t,p,rad]/n
                    endif
                 endif

                 if a[t,p,rad] eq 0 and t eq ntheta-1 and p eq nphi-1 then begin ; interpolating for theta = 357 deg and phi = 357 deg (Density)
                    for q=-1,0 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if u gt 0 and a[t+q,0,rad+g] ne 0 and a[t+q,0,rad+g] eq a[t+q,0,rad+g] then begin
                                a[t,p,rad]=a[t+q,0,rad+g]+a[t,p,rad]
                                na=na+1                             
                             endif else if u le 0 then begin
                                if a[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and a[t+q,p+u,rad+g] eq a[t+q,p+u,rad+g] then begin
                                   a[t,p,rad]=a[t+q,p+u,rad+g]+a[t,p,rad]
                                   na=na+1.
                                endif
                             endif
                          endfor
                       endfor
                    endfor
                    if na ne 0 then begin
                       a[t,p,rad]=a[t,p,rad]/na
                    endif
                 endif

                 if np_stream[t,p,rad] eq 0 and t eq ntheta-1 and p eq nphi-1 then begin ; interpolating for theta = 357 deg and phi = 357 deg (Density)
                    for q=-1,0 do begin
                       for u=-1,0 do begin
                          for g=-1,1 do begin
                             if np_stream[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and np_stream[t+q,p+u,rad+g] eq np_stream[t+q,p+u,rad+g] then begin
                                np_stream[t,p,rad]=np_stream[t+q,p+u,rad+g]+np_stream[t,p,rad]
                                tp_stream[t,p,rad]=tp_stream[t+q,p+u,rad+g]+tp_stream[t,p,rad]
                                vr_stream[t,p,rad]=vr_stream[t+q,p+u,rad+g]+vr_stream[t,p,rad]
                                vt_stream[t,p,rad]=vt_stream[t+q,p+u,rad+g]+vt_stream[t,p,rad]
                                nd=nd+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       np_stream[t,p,rad]=np_stream[t,p,rad]/nd
                       tp_stream[t,p,rad]=tp_stream[t,p,rad]/nd
                       vr_stream[t,p,rad]=vr_stream[t,p,rad]/nd
                       vt_stream[t,p,rad]=vt_stream[t,p,rad]/nd
                    endif
                 endif
                 
                 
                 if cool[t,p,rad] eq 0 and t eq ntheta-1 and p ne nphi-1 and p ne 0 then begin ; interpolating for theta = 357 deg and 3 deg < phi < 357 deg (Extinction)
                    for q=-1,0 do begin
                       for u=-1,1 do begin
                          if cool[t+q,p+u,rad] ne 0 and abs(q)+abs(u) ne 0 and cool[t+q,p+u,rad] eq cool[t+q,p+u,rad] then begin
                             cool[t,p,rad]=cool[t+q,p+u,rad]+cool[t,p,rad]
                             n=n+1.
                          endif
                       endfor
                    endfor
                    if n ne 0 then begin
                       cool[t,p,rad]=cool[t,p,rad]/n
                    endif
                 endif

                  if a[t,p,rad] eq 0 and t eq ntheta-1 and p ne nphi-1 and p ne 0 then begin ; interpolating for theta = 357 deg and 3 deg < phi < 357 deg (Density)
                    for q=-1,0 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if a[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and a[t+q,p+u,rad+g] eq a[t+q,p+u,rad+g] then begin
                                a[t,p,rad]=a[t+q,p+u,rad+g]+a[t,p,rad]
                                na=na+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if na ne 0 then begin
                       a[t,p,rad]=a[t,p,rad]/na
                    endif                  
                 endif

                  if np_stream[t,p,rad] eq 0 and t eq ntheta-1 and p ne nphi-1 and p ne 0 then begin ; interpolating for theta = 357 deg and 3 deg < phi < 357 deg (Density)
                    for q=-1,0 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if np_stream[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and np_stream[t+q,p+u,rad+g] eq np_stream[t+q,p+u,rad+g] then begin
                                np_stream[t,p,rad]=np_stream[t+q,p+u,rad+g]+np_stream[t,p,rad]
                                tp_stream[t,p,rad]=tp_stream[t+q,p+u,rad+g]+tp_stream[t,p,rad]
                                vr_stream[t,p,rad]=vr_stream[t+q,p+u,rad+g]+vr_stream[t,p,rad]
                                vt_stream[t,p,rad]=vt_stream[t+q,p+u,rad+g]+vt_stream[t,p,rad]
                                nd=nd+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       np_stream[t,p,rad]=np_stream[t,p,rad]/nd
                       tp_stream[t,p,rad]=tp_stream[t,p,rad]/nd
                       vr_stream[t,p,rad]=vr_stream[t,p,rad]/nd
                       vt_stream[t,p,rad]=vt_stream[t,p,rad]/nd
                    endif                  
                 endif
                  
                 if cool[t,p,rad] eq 0 and t ne ntheta-1 and p eq nphi-1 and t ne 0 then begin ; interpolating for phi = 357 deg and 3 deg < theta < 357 deg (Extinction)
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                          if u gt 0 and cool[t+q,0,rad] ne 0 and cool[t+q,0,rad] eq cool[t+q,0,rad] then begin
                             cool[t,p,rad]=cool[t+q,0,rad]+cool[t,p,rad]
                             n=n+1
                          endif else if u le 0 then begin
                             if cool[t+q,p+u,rad] ne 0 and abs(q)+abs(u) ne 0 and cool[t+q,p+u,rad] eq cool[t+q,p+u,rad] then begin
                                cool[t,p,rad]=cool[t+q,p+u,rad]+cool[t,p,rad]
                                n=n+1.
                             endif
                          endif
                        endfor
                    endfor
                    if n ne 0 then begin
                       cool[t,p,rad]=cool[t,p,rad]/n
                    endif
                 endif

                 if a[t,p,rad] eq 0 and t ne ntheta-1 and p eq nphi-1 and t ne 0 then begin ; interpolating for phi = 357 deg and 3 deg < theta < 357 deg (Density)
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if u gt 0 and a[t+q,0,rad+g] ne 0 and a[t+q,0,rad+g] eq a[t+q,0,rad+g] then begin
                                a[t,p,rad]=a[t+q,0,rad+g]+a[t,p,rad]
                                na=na+1
                                ;qp=1.
                                ;print, t+q, 0, rad+g, a[t+q,0,rad+g],ptra_ct[t,p,rad], '     new'
                             endif else if u le 0 then begin
                                if a[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and a[t+q,p+u,rad+g] eq a[t+q,p+u,rad+g] then begin
                                   a[t,p,rad]=a[t+q,p+u,rad+g]+a[t,p,rad]
                                   na=na+1.
                                   ;if qp eq 1. then print, t+q,p+u, rad+g, a[t+q,p+u,rad+g], ptra_ct[t,p,rad],'     old'
                                endif
                             endif
                          endfor
                       endfor
                    endfor
                    if na ne 0 then begin
                       a[t,p,rad]=a[t,p,rad]/na
                       ;if qp eq 1. then stop
                    endif 
                 endif

                 if np_stream[t,p,rad] eq 0 and t ne ntheta-1 and p eq nphi-1 and t ne 0 then begin ; interpolating for phi = 357 deg and 3 deg < theta < 357 deg (Density)
                    for q=-1,1 do begin
                       for u=-1,0 do begin
                          for g=-1,1 do begin
                             if np_stream[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and np_stream[t+q,p+u,rad+g] eq np_stream[t+q,p+u,rad+g] then begin
                                np_stream[t,p,rad]=np_stream[t+q,p+u,rad+g]+np_stream[t,p,rad]
                                tp_stream[t,p,rad]=tp_stream[t+q,p+u,rad+g]+tp_stream[t,p,rad]
                                vr_stream[t,p,rad]=vr_stream[t+q,p+u,rad+g]+vr_stream[t,p,rad]
                                vt_stream[t,p,rad]=vt_stream[t+q,p+u,rad+g]+vt_stream[t,p,rad]
                                nd=nd+1.
                              endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       np_stream[t,p,rad]=np_stream[t,p,rad]/nd
                       tp_stream[t,p,rad]=tp_stream[t,p,rad]/nd
                       vr_stream[t,p,rad]=vr_stream[t,p,rad]/nd
                       vt_stream[t,p,rad]=vt_stream[t,p,rad]/nd
                    endif 
                 endif
                
                 
                 if cool[t,p,rad] eq 0 and t eq ntheta-1 and p eq 0 then begin ; interpolating for theta = 357 deg and phi = 3 deg (Extinction)
                    for q=-1,0 do begin
                       for u=-1,1 do begin
                          if u lt 0 and cool[t+q,nphi-1,rad] ne 0 and cool[t+q,nphi-1,rad] eq cool[t+q,nphi-1,rad] then begin
                             cool[t,p,rad]=cool[t+q,nphi-1,rad]+cool[t,p,rad]
                             n=n+1
                          endif else if u ge 0 then begin
                             if cool[t+q,p+u,rad] ne 0 and abs(q)+abs(u) ne 0 and cool[t+q,p+u,rad] eq cool[t+q,p+u,rad] then begin
                                cool[t,p,rad]=cool[t+q,p+u,rad]+cool[t,p,rad]
                                n=n+1.
                             endif
                          endif
                       endfor
                    endfor
                    if n ne 0 then begin
                       cool[t,p,rad]=cool[t,p,rad]/n
                    endif
                 endif

                 if a[t,p,rad] eq 0 and t eq ntheta-1 and p eq 0 then begin ; interpolating for theta = 357 deg and phi = 3 deg (Density)
                    for q=-1,0 do begin
                       for u=-1,1 do begin             
                          for g=-1,1 do begin
                             if u lt 0 and a[t+q,nphi-1,rad+g] ne 0 and a[t+q,nphi-1,rad+g] eq a[t+q,nphi-1,rad+g] then begin
                                a[t,p,rad]=a[t+q,nphi-1,rad+g]+a[t,p,rad]
                                na=na+1  
                             endif else if u ge 0 then begin
                                if a[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and a[t+q,p+u,rad+g] eq a[t+q,p+u,rad+g] then begin
                                   a[t,p,rad]=a[t+q,p+u,rad+g]+a[t,p,rad]
                                   na=na+1.
                                endif
                             endif
                          endfor
                       endfor
                    endfor
                    if na ne 0 then begin
                       a[t,p,rad]=a[t,p,rad]/na
                    endif
                 endif

                 if np_stream[t,p,rad] eq 0 and t eq ntheta-1 and p eq 0 then begin ; interpolating for theta = 357 deg and phi = 3 deg (Density)
                    for q=-1,0 do begin
                       for u=0,1 do begin             
                          for g=-1,1 do begin
                             if np_stream[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and np_stream[t+q,p+u,rad+g] eq np_stream[t+q,p+u,rad+g] then begin
                                np_stream[t,p,rad]=np_stream[t+q,p+u,rad+g]+np_stream[t,p,rad]
                                tp_stream[t,p,rad]=tp_stream[t+q,p+u,rad+g]+tp_stream[t,p,rad]
                                vr_stream[t,p,rad]=vr_stream[t+q,p+u,rad+g]+vr_stream[t,p,rad]
                                vt_stream[t,p,rad]=vt_stream[t+q,p+u,rad+g]+vt_stream[t,p,rad]
                                nd=nd+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       np_stream[t,p,rad]=np_stream[t,p,rad]/nd
                       tp_stream[t,p,rad]=tp_stream[t,p,rad]/nd
                       vr_stream[t,p,rad]=vr_stream[t,p,rad]/nd
                       vt_stream[t,p,rad]=vt_stream[t,p,rad]/nd
                    endif
                 endif
                 
                 if cool[t,p,rad] eq 0 and t ne 0 and p ne 0 and t ne ntheta-1 and p ne nphi-1 then begin ; interpolating for 3 deg < theta < 357 deg and 3 deg < phi < 357 deg (Extinction)
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                             if cool[t+q,p+u,rad] ne 0 and abs(q)+abs(u) ne 0 and cool[t+q,p+u,rad] eq cool[t+q,p+u,rad] then begin
                                cool[t,p,rad]=cool[t+q,p+u,rad]+cool[t,p,rad]
                                n=n+1.
                             endif
                       endfor
                    endfor
                    if n ne 0 then begin
                       cool[t,p,rad]=cool[t,p,rad]/n
                    endif
                 endif
                 
                 if a[t,p,rad] eq 0 and t ne 0 and p ne 0 and t ne ntheta-1 and p ne nphi-1 then begin ; interpolating for 3 deg < theta < 357 deg and 3 deg < phi < 357 deg (Density)
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if a[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and a[t+q,p+u,rad+g] eq a[t+q,p+u,rad+g] then begin
                                a[t,p,rad]=a[t+q,p+u,rad+g]+a[t,p,rad]
                                na=na+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if na ne 0 then begin
                       a[t,p,rad]=a[t,p,rad]/na
                    endif
                 endif

                 if np_stream[t,p,rad] eq 0 and t ne 0 and p ne 0 and t ne ntheta-1 and p ne nphi-1 then begin ; interpolating for 3 deg < theta < 357 deg and 3 deg < phi < 357 deg (Density)
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if np_stream[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and np_stream[t+q,p+u,rad+g] eq np_stream[t+q,p+u,rad+g] then begin
                                np_stream[t,p,rad]=np_stream[t+q,p+u,rad+g]+np_stream[t,p,rad]
                                tp_stream[t,p,rad]=tp_stream[t+q,p+u,rad+g]+tp_stream[t,p,rad]
                                vr_stream[t,p,rad]=vr_stream[t+q,p+u,rad+g]+vr_stream[t,p,rad]
                                vt_stream[t,p,rad]=vt_stream[t+q,p+u,rad+g]+vt_stream[t,p,rad]
                                nd=nd+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       np_stream[t,p,rad]=np_stream[t,p,rad]/nd
                       tp_stream[t,p,rad]=tp_stream[t,p,rad]/nd
                       vr_stream[t,p,rad]=vr_stream[t,p,rad]/nd
                       vt_stream[t,p,rad]=vt_stream[t,p,rad]/nd
                    endif
                 endif

              endif else cool[t,p,rad]=0.
           endfor
        endif    
     endfor
  endfor
endfor  

for i=0,1 do begin
   for rad=0,rupper do begin
      if i eq 0 and boundaries_arr[0,0,rad] eq 3 then begin
         w=where(boundaries_arr[1,*,rad] eq 3, ct)
         if w[0] ge 0 then cool[0,*,rad]=mean(cool[1,w,rad])
         if w[0] ge 0 then a[0,*,rad]=mean(a[1,w,rad])
         if w[0] ge 0 then np_stream[0,*,rad]=mean(np_stream[1,w,rad])
         if w[0] ge 0 then tp_stream[0,*,rad]=mean(tp_stream[1,w,rad])
         if w[0] ge 0 then vr_stream[0,*,rad]=mean(vr_stream[1,w,rad])
         if w[0] ge 0 then vt_stream[0,*,rad]=mean(vt_stream[1,w,rad])
      endif else if i eq 1 and boundaries_arr[0,0,rad] eq 3 then begin
         if cool[0,0,rad] eq 0 and cool[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then cool[0,*,rad]=cool[0,0,rad+1]
         if cool[0,0,rad] eq 0 and cool[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then cool[0,*,rad]=cool[0,0,rad-1]
         if a[0,0,rad] eq 0 and a[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then a[0,0,rad]=a[0,*,rad+1]
         if a[0,0,rad] eq 0 and a[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then a[0,0,rad]=a[0,*,rad-1]
         if np_stream[0,0,rad] eq 0 and np_stream[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then np_stream[0,0,rad]=np_stream[0,*,rad+1]
         if np_stream[0,0,rad] eq 0 and np_stream[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then np_stream[0,0,rad]=np_stream[0,*,rad-1]
         if tp_stream[0,0,rad] eq 0 and tp_stream[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then tp_stream[0,0,rad]=tp_stream[0,*,rad+1]
         if tp_stream[0,0,rad] eq 0 and tp_stream[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then tp_stream[0,0,rad]=tp_stream[0,*,rad-1]
         if vr_stream[0,0,rad] eq 0 and vr_stream[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then vr_stream[0,0,rad]=vr_stream[0,*,rad+1]
         if vr_stream[0,0,rad] eq 0 and vr_stream[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then vr_stream[0,0,rad]=vr_stream[0,*,rad-1]
         if vt_stream[0,0,rad] eq 0 and vt_stream[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then vt_stream[0,0,rad]=vt_stream[0,*,rad+1]
         if vt_stream[0,0,rad] eq 0 and vt_stream[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then vt_stream[0,0,rad]=vt_stream[0,*,rad-1]
        
      endif
   endfor
endfor
  
  vx[0,0,*]=vx[0,1,*]
  vy[0,0,*]=vy[0,1,*]
  vz[0,0,*]=vz[0,1,*]


end
