pro read_stream_g, lc, cool, vx, vy, vz, ntr, nref, nen, etr, eref, een, em, kappa_ref, vshock, esw
  
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
  ntr0=dblarr(num-1)
  nref0=dblarr(num-1)
  nen0=dblarr(num-1)
  etr0=dblarr(num-1)
  eref0=dblarr(num-1)
  een0=dblarr(num-1)
  em0=dblarr(num-1)
  kappa_ref0=dblarr(num-1)
  vshock0=dblarr(num-1)
  esw0=dblarr(num-1)
  ptr_ct=dblarr(ntheta,nphi,nr)
  ptrn_ct=dblarr(ntheta,nphi,nr)
  btest=dblarr(ntheta,nphi,nr)

;fixing phi
  phi=atan(y,x)*57.2958
  for q=0, num-1 do if phi[q] lt 0 then phi[q]=360.-abs(phi[q])
  
  v=sqrt(vr^2.+vt^2.)
  cs=sqrt(5*(2*tp)*kb/(3*mp)) ; multiply by 2 to account for electrons
  M2=(v^2.)/(cs^2.)

  dr0=r[1:num-1]-r[0:num-2]
  dr0=dr0*au_cm ; convert from au to cm
  dtheta0=theta[1:num-1]-theta[0:num-2]
  dtheta0=dtheta0*(!pi/180.) ; convert from deg to rad
  dphi0=phi[1:num-1]-phi[0:num-2]
  dphi0=dtheta0*(!pi/180.) ; convert from deg to rad
  ds=sqrt(dr0^2.+r[1:num-1]^2.*dtheta0^2.+r[1:num-1]^2.*sin(theta[1:num-1]*(!pi/180.))^2.*dphi0^2.) ; spherical polar line element

 tau=-(nh[0:num-2]/v[1:num-1])*ds[0:num-2]

if moscow_regs eq 1 then begin
  for i=1,num-2 do begin
     if reg[i-1] eq 2 and reg[i] lt 2 and reg[i] gt 1 then reg[i]=2.
     if reg[i] ge 1 and reg[i] lt 2 and reg[i+1] eq 2  and r[i] lt 200. then begin
        ; 4 population case
        r_ts0=r[i]/82. ; normalized to V2
        if r_ts0 ge 1. then r_ts=alog(r_ts0) else r_ts=alog(1.); putting into log space and fixing for r_ts less than V2 (limiting by hybrid)
        if r_ts0 gt 1.62 then r_ts=alog(1.62) ; fixing for tail points (limited by hybrid ratio of tail-to-V2)
        ntr0[i]=exp(-4.65496*r_ts^2+2.46160*r_ts-1.71480)
        nref0[i]=exp(-2.27104*r_ts^2+0.395963*r_ts-2.65926)
        nen0[i]=exp(8.10281*r_ts^2-3.06275*r_ts-2.52573)
        etr0[i]=exp(0.101026*r_ts^2+0.472896*r_ts-0.579818)
        eref0[i]=exp(-15.8367*r_ts^2+3.91212*r_ts-1.71480)
        een0[i]=exp(2.73671*r_ts^2-0.477155*r_ts-0.223144)
        em0[i]=exp(18.7905*r_ts^2-10.5773*r_ts+3.21888)
        kappa_ref0[i]=exp(-13.5408*r_ts^2+6.52448*r_ts+2.01490); - this is for delta=7.5 @ V2
        esw0[i]=exp(2.73671*r_ts^2-0.477155*r_ts-3.21888)
        ;if r[i]/82. lt 1. then nen0[i]=1-0.75-ntr0[i]-nref0[i]; Full: 1-0.70-ntr0[i]-nref0[i]
        ;if r[i]/82. lt 1. then een0[i]=1-0.04-etr0[i]-eref0[i]; Full: 1-0.06-etr0[i]-eref0[i]
        ;if 1-ntr0[i]-nref0[i] le 0.65 then begin
	;	ndiff=1-ntr0[i]-nref0[i]
        ;        nswdiff=0.65-ndiff
	;	ntr0[i]=ntr0[i]-nswdiff/2.
	;	nref0[i]=nref0[i]-nswdiff/2.
	;endif
	;if 1-etr0[i]-eref0[i] le 0.04 then begin
	;	ediff=1-etr0[i]-eref0[i]
        ;        eswdiff=0.04-ediff
        ;        etr0[i]=etr0[i]-eswdiff/2.
        ;        eref0[i]=eref0[i]-eswdiff/2.
	;endif
    endif else if reg[i] eq 2 and r[i] eq r[i] and r[i-1] eq r[i-1] then begin
        ntr0[i]=ntr0[i-1]
        nref0[i]=nref0[i-1]
	nen0[i]=nen0[i-1]
        etr0[i]=etr0[i-1]
        eref0[i]=eref0[i-1]
	een0[i]=een0[i-1]
        em0[i]=em0[i-1]
        kappa_ref0[i]=kappa_ref0[i-1]
        esw0[i]=esw0[i-1]
        if ntr0[i] gt 1 or nref0[i] gt 1 or etr0[i] gt 1 or eref0[i] gt 1 then stop
    endif
    if reg[i] eq 2 and r[i] eq r[i] and r[i-1] eq r[i-1] and vshock0[i-1] eq 0. then vshock0[i]=v[i]
    if reg[i] eq 2 and r[i] eq r[i] and r[i-1] eq r[i-1] and vshock0[i-1] gt 0. then vshock0[i]=vshock0[i-1]
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
     if M2[i-1] gt 0.64 and M2[i] lt 0.64 and r[i] lt 130 then begin
       ; 4 population case
        r_ts0=r[i]/78.5;69.;78.5 ; normalized to V2
        if r_ts0 ge 1. then r_ts=alog(r_ts0) else r_ts=alog(1.); putting into log space and fixing for r_ts less than V2 (limiting by hybrid)
        if r_ts0 gt 1.62 then r_ts=alog(1.62) ; fixing for tail points (limited by hybrid ratio of tail-to-V2) 
        ntr0[i]=exp(-4.65496*r_ts^2+2.46160*r_ts-1.71480)
        nref0[i]=exp(-2.27104*r_ts^2+0.395963*r_ts-2.65926)
        nen0[i]=exp(8.10281*r_ts^2-3.06275*r_ts-2.52573)
        etr0[i]=exp(0.101026*r_ts^2+0.472896*r_ts-0.579818)
        eref0[i]=exp(-15.8367*r_ts^2+3.91212*r_ts-1.71480)
        een0[i]=exp(2.73671*r_ts^2-0.477155*r_ts-0.223144)
        em0[i]=exp(18.7905*r_ts^2-10.5773*r_ts+3.21888)
        kappa_ref0[i]=exp(-13.5408*r_ts^2+6.52448*r_ts+2.01490); - this is for delta=7.5 @ V2
        esw0[i]=exp(2.73671*r_ts^2-0.477155*r_ts-3.21888)
        ;if r[i]/78.5 lt 1. then nen0[i]=1-0.75-ntr0[i]-nref0[i]; Full: 1-0.70-ntr0[i]-nref0[i]
        ;if r[i]/78.5 lt 1. then een0[i]=1-0.04-etr0[i]-eref0[i]; Full: 1-0.06-etr0[i]-eref0[i]
        ;if 1-ntr0[i]-nref0[i] le 0.65 then begin
        ;        ndiff=1-ntr0[i]-nref0[i]
        ;        nswdiff=0.65-ndiff
        ;        ntr0[i]=ntr0[i]-nswdiff/2.
        ;        nref0[i]=nref0[i]-nswdiff/2.
        ;endif
        ;if 1-etr0[i]-eref0[i] le 0.04 then begin
        ;        ediff=1-etr0[i]-eref0[i]
        ;        eswdiff=0.04-ediff
        ;        etr0[i]=etr0[i]-eswdiff/2.
        ;        eref0[i]=eref0[i]-eswdiff/2.
        ;endif
     endif else if M2[i] lt 0.64 and r[i] eq r[i] and r[i-1] eq r[i-1] then begin
       ntr0[i]=ntr0[i-1]
       nref0[i]=nref0[i-1]
       nen0[i]=nen0[i-1]
       etr0[i]=etr0[i-1]
       eref0[i]=eref0[i-1]
       een0[i]=een0[i-1]
       em0[i]=em0[i-1]
       kappa_ref0[i]=kappa_ref0[i-1]
       esw0[i]=esw0[i-1]
       if ntr0[i] gt 1 or nref0[i] gt 1 or etr0[i] gt 1 or eref0[i] gt 1 then stop
     endif
     if M2[i] lt 0.64 and r[i] eq r[i] and r[i-1] eq r[i-1] and vshock0[i-1] eq 0. then vshock0[i]=v[i]
     if M2[i] lt 0.64 and r[i] eq r[i] and r[i-1] eq r[i-1] and vshock0[i-1] gt 0. then vshock0[i]=vshock0[i-1]
  endfor  
;stop
  nan=0./0.
; Ensuring only considering points starting from TS
  for i=1, num-2 do begin
     if M2[i] gt 0.81 then tau[i]=0
     if (r[i] eq r[i]) and (r[i-1] eq r[i-1]) then ext[i]=tau[i]+ext[i-1] $
     else ext[i]=0
  endfor
endelse

  if num lt 2000 then streamline_test, ext, r, nh, v, vr, vt, np, tp

  dist=where(ext gt -1.01 and ext lt -0.99,ct)
  r_in=round((r-r_lower)/dr) ; indices for radius
  t_in=round((theta-2.9)/dtheta) ; indices for theta
  p_in=round((phi-2.9)/dphi)     ; indices for phi

 rlower=(r_lower-r_inner)/dr; 23
 rupper=(r_upper-r_inner)/dr

for i=0,10 do begin
  for t=0,ntheta-1 do begin
     for p=0,nphi-1 do begin
        if i eq 0 then begin ; putting extinction values into 3D array of theta, phi and r
           print, string((nphi*t+p)/(nphi*ntheta)*100.)+'%'
           for rad=rlower, rupper do begin
             if boundaries_arr[t,p,rad] eq 3. or boundaries_arr[t,p,rad] eq 6 then begin
                 ptr0=where(p_in eq p and t_in eq t and r_in eq rad,ptr0ct) ; Actually 'ptr-1' based on the way tau and ext is calculated
		 if ptr0ct gt 0 then ptr_test=where(ext[ptr0] eq ext[ptr0],ptrct) else ptrct=0
                 if ptrct gt 0 then ptr=ptr0[ptr_test] else ptr=-1
		 if ptrct gt 0 then ptr_ct[t,p,rad]=ptrct
                 if ptrct gt 0 then cool[t,p,rad]=mean(ext[ptr])
                 if cool[t,p,rad] ne cool[t,p,rad] and t eq 0 and p eq 0 then cool[t,p,rad]=0.
                 if ptr(0) ge 0 then vx[t,p,rad]=mean(ux[ptr]) else vx[t,p,rad]=0. ; for printing out streamlines
                 if ptr(0) ge 0 then vy[t,p,rad]=mean(uy[ptr]) else vy[t,p,rad]=0.
                 if ptr(0) ge 0 then vz[t,p,rad]=mean(uz[ptr]) else vz[t,p,rad]=0.

                 if ptr0ct gt 0 then ptrn_test=where(ntr0[ptr0] eq ntr0[ptr0] and ntr0[ptr0] ne 0,ptrnct) else ptrnct=0
                 if ptrnct gt 0 then ptrn=ptr0[ptrn_test] else ptrn=-1
		 if ptrnct gt 0 then ptrn_ct[t,p,rad]=ptrnct
                 if ptrn(0) ge 0 then ntr[t,p,rad]=mean(ntr0[ptrn])
                 if ptrn(0) ge 0 then nref[t,p,rad]=mean(nref0[ptrn])
                 if ptrn(0) ge 0 then etr[t,p,rad]=mean(etr0[ptrn])
                 if ptrn(0) ge 0 then eref[t,p,rad]=mean(eref0[ptrn])
                 if ptrn(0) ge 0 then nen[t,p,rad]=mean(nen0[ptrn])
                 if ptrn(0) ge 0 then een[t,p,rad]=mean(een0[ptrn])
                 if ptrn(0) ge 0 then em[t,p,rad]=mean(em0[ptrn])
                 if ptrn(0) ge 0 then kappa_ref[t,p,rad]=mean(kappa_ref0[ptrn])
                 if ptrn(0) ge 0 then vshock[t,p,rad]=mean(vshock0[ptrn])
                 if ptrn(0) ge 0 then esw[t,p,rad]=mean(esw0[ptrn])
              endif else cool[t,p,rad]=0.                 
           endfor
        endif

        if i eq 1 then begin    ; first wave of interpolating zero extinction values
           for rad=rlower, rupper do begin
              if boundaries_arr[t,p,rad] eq 3. or boundaries_arr[t,p,rad] eq 6 then begin
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

                 if ntr[t,p,rad] eq 0 and t ne 0 and p ne 0 and t ne ntheta-1 and p ne nphi-1 then begin ; interpolating for 3 deg < theta < 357 deg and 3 deg < phi < 357 deg
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if ntr[t+q,p+u,rad+g] ne 0 and ptrn_ct[t+q,p+u,rad+g] gt 0 and abs(q)+abs(u)+abs(g) ne 0 and ntr[t+q,p+u,rad+g] eq ntr[t+q,p+u,rad+g] then begin
                                ntr[t,p,rad]=ntr[t+q,p+u,rad+g]+ntr[t,p,rad]
                                nref[t,p,rad]=nref[t+q,p+u,rad+g]+nref[t,p,rad]
                                etr[t,p,rad]=etr[t+q,p+u,rad+g]+etr[t,p,rad]
                                eref[t,p,rad]=eref[t+q,p+u,rad+g]+eref[t,p,rad]
                                nen[t,p,rad]=nen[t+q,p+u,rad+g]+nen[t,p,rad]
                                een[t,p,rad]=een[t+q,p+u,rad+g]+een[t,p,rad]
                                em[t,p,rad]=em[t+q,p+u,rad+g]+em[t,p,rad]
                                kappa_ref[t,p,rad]=kappa_ref[t+q,p+u,rad+g]+kappa_ref[t,p,rad]
                                vshock[t,p,rad]=vshock[t+q,p+u,rad+g]+vshock[t,p,rad]
                                esw[t,p,rad]=esw[t+q,p+u,rad+g]+esw[t,p,rad] 
                               nd=nd+1. ; if nonzero extinction, adds to index for averaging
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
		       ntr[t,p,rad]=ntr[t,p,rad]/nd
		       nref[t,p,rad]=nref[t,p,rad]/nd
		       etr[t,p,rad]=etr[t,p,rad]/nd
		       eref[t,p,rad]=eref[t,p,rad]/nd
                       nen[t,p,rad]=nen[t,p,rad]/nd
                       een[t,p,rad]=een[t,p,rad]/nd
                       em[t,p,rad]=em[t,p,rad]/nd
                       kappa_ref[t,p,rad]=kappa_ref[t,p,rad]/nd
                       vshock[t,p,rad]=vshock[t,p,rad]/nd
                    endif
                 endif
              endif else cool[t,p,rad]=0. 
           endfor
        endif  
        
        if i ge 2 then begin
           for rad=rlower, rupper do begin
              if boundaries_arr[t,p,rad] eq 3. or boundaries_arr[t,p,rad] eq 6. then begin
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

                 if ntr[t,p,rad] eq 0 and t ne 0 and p eq 0 and t ne ntheta-1 then begin ; interpolating for phi = 3 deg and 3 deg < theta < 357 deg (Density)
                    for q=-1,1 do begin
                       for u=0,1 do begin
                          for g=-1,1 do begin
                             if ntr[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and ntr[t+q,p+u,rad+g] eq ntr[t+q,p+u,rad+g] then begin
                                ntr[t,p,rad]=ntr[t+q,p+u,rad+g]+ntr[t,p,rad]
                                nref[t,p,rad]=nref[t+q,p+u,rad+g]+nref[t,p,rad]
                                etr[t,p,rad]=etr[t+q,p+u,rad+g]+etr[t,p,rad]
                                eref[t,p,rad]=eref[t+q,p+u,rad+g]+eref[t,p,rad]
                                nen[t,p,rad]=nen[t+q,p+u,rad+g]+nen[t,p,rad]
                                een[t,p,rad]=een[t+q,p+u,rad+g]+een[t,p,rad]
                                em[t,p,rad]=em[t+q,p+u,rad+g]+em[t,p,rad]
                                kappa_ref[t,p,rad]=kappa_ref[t+q,p+u,rad+g]+kappa_ref[t,p,rad]
                                vshock[t,p,rad]=vshock[t+q,p+u,rad+g]+vshock[t,p,rad]
                                esw[t,p,rad]=esw[t+q,p+u,rad+g]+esw[t,p,rad]
                                nd=nd+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       ntr[t,p,rad]=ntr[t,p,rad]/nd
                       nref[t,p,rad]=nref[t,p,rad]/nd
                       etr[t,p,rad]=etr[t,p,rad]/nd
                       eref[t,p,rad]=eref[t,p,rad]/nd
                       nen[t,p,rad]=nen[t,p,rad]/nd
                       een[t,p,rad]=een[t,p,rad]/nd
                       em[t,p,rad]=em[t,p,rad]/nd
                       kappa_ref[t,p,rad]=kappa_ref[t,p,rad]/nd
                       vshock[t,p,rad]=vshock[t,p,rad]/nd
                       esw[t,p,rad]=esw[t,p,rad]/nd
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

                 if ntr[t,p,rad] eq 0 and t eq ntheta-1 and p eq nphi-1 then begin ; interpolating for theta = 357 deg and phi = 357 deg (Density)
                    for q=-1,0 do begin
                       for u=-1,0 do begin
                          for g=-1,1 do begin
                             if ntr[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and ntr[t+q,p+u,rad+g] eq ntr[t+q,p+u,rad+g] then begin
                                ntr[t,p,rad]=ntr[t+q,p+u,rad+g]+ntr[t,p,rad]
                                nref[t,p,rad]=nref[t+q,p+u,rad+g]+nref[t,p,rad]
                                etr[t,p,rad]=etr[t+q,p+u,rad+g]+etr[t,p,rad]
                                eref[t,p,rad]=eref[t+q,p+u,rad+g]+eref[t,p,rad]
                                nen[t,p,rad]=nen[t+q,p+u,rad+g]+nen[t,p,rad]
                                een[t,p,rad]=een[t+q,p+u,rad+g]+een[t,p,rad]
                                em[t,p,rad]=em[t+q,p+u,rad+g]+em[t,p,rad]
                                kappa_ref[t,p,rad]=kappa_ref[t+q,p+u,rad+g]+kappa_ref[t,p,rad]
                                vshock[t,p,rad]=vshock[t+q,p+u,rad+g]+vshock[t,p,rad]
                                esw[t,p,rad]=esw[t+q,p+u,rad+g]+esw[t,p,rad]
                                nd=nd+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       ntr[t,p,rad]=ntr[t,p,rad]/nd
                       nref[t,p,rad]=nref[t,p,rad]/nd
                       etr[t,p,rad]=etr[t,p,rad]/nd
                       eref[t,p,rad]=eref[t,p,rad]/nd
                       nen[t,p,rad]=nen[t,p,rad]/nd
                       een[t,p,rad]=een[t,p,rad]/nd
                       em[t,p,rad]=em[t,p,rad]/nd
                       kappa_ref[t,p,rad]=kappa_ref[t,p,rad]/nd
                       vshock[t,p,rad]=vshock[t,p,rad]/nd
                       esw[t,p,rad]=esw[t,p,rad]/nd
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

                  if ntr[t,p,rad] eq 0 and t eq ntheta-1 and p ne nphi-1 and p ne 0 then begin ; interpolating for theta = 357 deg and 3 deg < phi < 357 deg (Density)
                    for q=-1,0 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if ntr[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and ntr[t+q,p+u,rad+g] eq ntr[t+q,p+u,rad+g] then begin
                                ntr[t,p,rad]=ntr[t+q,p+u,rad+g]+ntr[t,p,rad]
                                nref[t,p,rad]=nref[t+q,p+u,rad+g]+nref[t,p,rad]
                                etr[t,p,rad]=etr[t+q,p+u,rad+g]+etr[t,p,rad]
                                eref[t,p,rad]=eref[t+q,p+u,rad+g]+eref[t,p,rad]
                                nen[t,p,rad]=nen[t+q,p+u,rad+g]+nen[t,p,rad]
                                een[t,p,rad]=een[t+q,p+u,rad+g]+een[t,p,rad]
                                em[t,p,rad]=em[t+q,p+u,rad+g]+em[t,p,rad]
                                kappa_ref[t,p,rad]=kappa_ref[t+q,p+u,rad+g]+kappa_ref[t,p,rad]
                                vshock[t,p,rad]=vshock[t+q,p+u,rad+g]+vshock[t,p,rad]
                                esw[t,p,rad]=esw[t+q,p+u,rad+g]+esw[t,p,rad]
                                nd=nd+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       ntr[t,p,rad]=ntr[t,p,rad]/nd
                       nref[t,p,rad]=nref[t,p,rad]/nd
                       etr[t,p,rad]=etr[t,p,rad]/nd
                       eref[t,p,rad]=eref[t,p,rad]/nd
                       nen[t,p,rad]=nen[t,p,rad]/nd
                       een[t,p,rad]=een[t,p,rad]/nd
                       em[t,p,rad]=em[t,p,rad]/nd
                       kappa_ref[t,p,rad]=kappa_ref[t,p,rad]/nd
                       vshock[t,p,rad]=vshock[t,p,rad]/nd
                       esw[t,p,rad]=esw[t,p,rad]/nd
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

                 if ntr[t,p,rad] eq 0 and t ne ntheta-1 and p eq nphi-1 and t ne 0 then begin ; interpolating for phi = 357 deg and 3 deg < theta < 357 deg (Density)
                    for q=-1,1 do begin
                       for u=-1,0 do begin
                          for g=-1,1 do begin
                             if ntr[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and ntr[t+q,p+u,rad+g] eq ntr[t+q,p+u,rad+g] then begin
                                ntr[t,p,rad]=ntr[t+q,p+u,rad+g]+ntr[t,p,rad]
                                nref[t,p,rad]=nref[t+q,p+u,rad+g]+nref[t,p,rad]
                                etr[t,p,rad]=etr[t+q,p+u,rad+g]+etr[t,p,rad]
                                eref[t,p,rad]=eref[t+q,p+u,rad+g]+eref[t,p,rad]
                                nen[t,p,rad]=nen[t+q,p+u,rad+g]+nen[t,p,rad]
                                een[t,p,rad]=een[t+q,p+u,rad+g]+een[t,p,rad]
                                em[t,p,rad]=em[t+q,p+u,rad+g]+em[t,p,rad]
                                kappa_ref[t,p,rad]=kappa_ref[t+q,p+u,rad+g]+kappa_ref[t,p,rad]
                                vshock[t,p,rad]=vshock[t+q,p+u,rad+g]+vshock[t,p,rad]
                                esw[t,p,rad]=esw[t+q,p+u,rad+g]+esw[t,p,rad]
                                nd=nd+1.
                              endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       ntr[t,p,rad]=ntr[t,p,rad]/nd
                       nref[t,p,rad]=nref[t,p,rad]/nd
                       etr[t,p,rad]=etr[t,p,rad]/nd
                       eref[t,p,rad]=eref[t,p,rad]/nd
                       nen[t,p,rad]=nen[t,p,rad]/nd
                       een[t,p,rad]=een[t,p,rad]/nd
                       em[t,p,rad]=em[t,p,rad]/nd
                       kappa_ref[t,p,rad]=kappa_ref[t,p,rad]/nd
                       vshock[t,p,rad]=vshock[t,p,rad]/nd
                       esw[t,p,rad]=esw[t,p,rad]/nd
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

                 if ntr[t,p,rad] eq 0 and t eq ntheta-1 and p eq 0 then begin ; interpolating for theta = 357 deg and phi = 3 deg (Density)
                    for q=-1,0 do begin
                       for u=0,1 do begin             
                          for g=-1,1 do begin
                             if ntr[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and ntr[t+q,p+u,rad+g] eq ntr[t+q,p+u,rad+g] then begin
                                ntr[t,p,rad]=ntr[t+q,p+u,rad+g]+ntr[t,p,rad]
                                nref[t,p,rad]=nref[t+q,p+u,rad+g]+nref[t,p,rad]
                                etr[t,p,rad]=etr[t+q,p+u,rad+g]+etr[t,p,rad]
                                eref[t,p,rad]=eref[t+q,p+u,rad+g]+eref[t,p,rad]
                                nen[t,p,rad]=nen[t+q,p+u,rad+g]+nen[t,p,rad]
                                een[t,p,rad]=een[t+q,p+u,rad+g]+een[t,p,rad]
                                em[t,p,rad]=em[t+q,p+u,rad+g]+em[t,p,rad]
                                kappa_ref[t,p,rad]=kappa_ref[t+q,p+u,rad+g]+kappa_ref[t,p,rad]
                                vshock[t,p,rad]=vshock[t+q,p+u,rad+g]+vshock[t,p,rad]
                                esw[t,p,rad]=esw[t+q,p+u,rad+g]+esw[t,p,rad]
                                nd=nd+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       ntr[t,p,rad]=ntr[t,p,rad]/nd
                       nref[t,p,rad]=nref[t,p,rad]/nd
                       etr[t,p,rad]=etr[t,p,rad]/nd
                       eref[t,p,rad]=eref[t,p,rad]/nd
                       nen[t,p,rad]=nen[t,p,rad]/nd
                       een[t,p,rad]=een[t,p,rad]/nd
                       em[t,p,rad]=em[t,p,rad]/nd
                       kappa_ref[t,p,rad]=kappa_ref[t,p,rad]/nd
                       vshock[t,p,rad]=vshock[t,p,rad]/nd
                       esw[t,p,rad]=esw[t,p,rad]/nd
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

                 if ntr[t,p,rad] eq 0 and t ne 0 and p ne 0 and t ne ntheta-1 and p ne nphi-1 then begin ; interpolating for 3 deg < theta < 357 deg and 3 deg < phi < 357 deg (Density)
                    for q=-1,1 do begin
                       for u=-1,1 do begin
                          for g=-1,1 do begin
                             if ntr[t+q,p+u,rad+g] ne 0 and abs(q)+abs(u)+abs(g) ne 0 and ntr[t+q,p+u,rad+g] eq ntr[t+q,p+u,rad+g] then begin
                                ntr[t,p,rad]=ntr[t+q,p+u,rad+g]+ntr[t,p,rad]
                                nref[t,p,rad]=nref[t+q,p+u,rad+g]+nref[t,p,rad]
                                etr[t,p,rad]=etr[t+q,p+u,rad+g]+etr[t,p,rad]
                                eref[t,p,rad]=eref[t+q,p+u,rad+g]+eref[t,p,rad]
                                nen[t,p,rad]=nen[t+q,p+u,rad+g]+nen[t,p,rad]
                                een[t,p,rad]=een[t+q,p+u,rad+g]+een[t,p,rad]
                                em[t,p,rad]=em[t+q,p+u,rad+g]+em[t,p,rad]
                                kappa_ref[t,p,rad]=kappa_ref[t+q,p+u,rad+g]+kappa_ref[t,p,rad]
                                vshock[t,p,rad]=vshock[t+q,p+u,rad+g]+vshock[t,p,rad] 
                                esw[t,p,rad]=esw[t+q,p+u,rad+g]+esw[t,p,rad]
                                nd=nd+1.
                             endif
                          endfor
                       endfor
                    endfor
                    if nd ne 0 then begin
                       ntr[t,p,rad]=ntr[t,p,rad]/nd
                       nref[t,p,rad]=nref[t,p,rad]/nd
                       etr[t,p,rad]=etr[t,p,rad]/nd
                       eref[t,p,rad]=eref[t,p,rad]/nd
                       nen[t,p,rad]=nen[t,p,rad]/nd
                       een[t,p,rad]=een[t,p,rad]/nd
                       em[t,p,rad]=em[t,p,rad]/nd
                       kappa_ref[t,p,rad]=kappa_ref[t,p,rad]/nd
                       vshock[t,p,rad]=vshock[t,p,rad]/nd
                       esw[t,p,rad]=esw[t,p,rad]/nd
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
         if w[0] ge 0 then ntr[0,*,rad]=mean(ntr[1,w,rad])
         if w[0] ge 0 then nref[0,*,rad]=mean(nref[1,w,rad])
         if w[0] ge 0 then etr[0,*,rad]=mean(etr[1,w,rad])
         if w[0] ge 0 then eref[0,*,rad]=mean(eref[1,w,rad])
         if w[0] ge 0 then nen[0,*,rad]=mean(nen[1,w,rad])
         if w[0] ge 0 then een[0,*,rad]=mean(een[1,w,rad])
         if w[0] ge 0 then em[0,*,rad]=mean(em[1,w,rad])
         if w[0] ge 0 then kappa_ref[0,*,rad]=mean(kappa_ref[1,w,rad])
         if w[0] ge 0 then vshock[0,*,rad]=mean(vshock[1,w,rad])
         if w[0] ge 0 then esw[0,*,rad]=mean(esw[1,w,rad])
      endif else if i eq 1 and boundaries_arr[0,0,rad] eq 3 then begin
         if cool[0,0,rad] eq 0 and cool[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then cool[0,*,rad]=cool[0,0,rad+1]
         if cool[0,0,rad] eq 0 and cool[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then cool[0,*,rad]=cool[0,0,rad-1]
         if ntr[0,0,rad] eq 0 and ntr[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then ntr[0,0,rad]=ntr[0,*,rad+1]
         if ntr[0,0,rad] eq 0 and ntr[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then ntr[0,0,rad]=ntr[0,*,rad-1]
         if nref[0,0,rad] eq 0 and nref[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then nref[0,0,rad]=nref[0,*,rad+1]
         if nref[0,0,rad] eq 0 and nref[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then nref[0,0,rad]=nref[0,*,rad-1]
         if etr[0,0,rad] eq 0 and etr[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then etr[0,0,rad]=etr[0,*,rad+1]
         if etr[0,0,rad] eq 0 and etr[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then etr[0,0,rad]=etr[0,*,rad-1]
         if eref[0,0,rad] eq 0 and eref[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then eref[0,0,rad]=eref[0,*,rad+1]
         if eref[0,0,rad] eq 0 and eref[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then eref[0,0,rad]=eref[0,*,rad-1]
         if nen[0,0,rad] eq 0 and nen[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then nen[0,0,rad]=nen[0,*,rad+1]
         if nen[0,0,rad] eq 0 and nen[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then nen[0,0,rad]=nen[0,*,rad-1]
         if een[0,0,rad] eq 0 and een[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then een[0,0,rad]=een[0,*,rad+1]
         if een[0,0,rad] eq 0 and een[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then een[0,0,rad]=een[0,*,rad-1]
         if em[0,0,rad] eq 0 and em[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then em[0,0,rad]=em[0,*,rad+1]
         if em[0,0,rad] eq 0 and em[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then em[0,0,rad]=em[0,*,rad-1]
         if kappa_ref[0,0,rad] eq 0 and kappa_ref[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then kappa_ref[0,0,rad]=kappa_ref[0,*,rad+1]
         if kappa_ref[0,0,rad] eq 0 and kappa_ref[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then kappa_ref[0,0,rad]=kappa_ref[0,*,rad-1]        
         if vshock[0,0,rad] eq 0 and vshock[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then vshock[0,0,rad]=vshock[0,*,rad+1]
         if vshock[0,0,rad] eq 0 and vshock[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then vshock[0,0,rad]=vshock[0,*,rad-1]
         if esw[0,0,rad] eq 0 and esw[0,0,rad+1] ne 0 and boundaries_arr[0,0,rad+1] eq 3 then esw[0,0,rad]=esw[0,*,rad+1]
         if esw[0,0,rad] eq 0 and esw[0,0,rad-1] ne 0 and boundaries_arr[0,0,rad-1] eq 3 then esw[0,0,rad]=esw[0,*,rad-1]
      endif
   endfor
endfor
  
  vx[0,0,*]=vx[0,1,*]
  vy[0,0,*]=vy[0,1,*]
  vz[0,0,*]=vz[0,1,*]

end
