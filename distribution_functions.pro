function kappa_fn,T,np,v,kappa0 ; kappa added by MZK
;Inputs:
; T = Temperature (Kelvin)
; np = number density of protons (cm^-3)
; v = velocity in reference frame of bulk speed 

common kappa_value, kappa
common constants_cgs,mp,kb,kev_erg,AU_cm
w=((2.*(kappa0-1.5)/kappa0)*kb*T/mp)^(1./2.) ;most probable speed in RF of bulk speed from pg. 44 of log book
f1=np/(w^3)/(!pi)^(3./2.)
f2=gamma(kappa0+1.)/gamma(kappa0-.5)/(kappa0)^(3./2.)
;f1=np/(2*!pi)/(kappa*w^2.)^(3./2.); added by MZK to agree with Pierrard & Lazar (2010)
;f2=gamma(kappa+1.)/gamma(kappa-.5)/gamma(3./2.); added by MZK to agree with Pierrard & Lazar (2010)
f3=(1+v^2./kappa0/(w^2))^(-1.-kappa0)
fk=f1*f2*f3

if fk ne fk then stop

return,fk
end
;------------------------------
function maxwellian_fn,T,np,v
;Inputs:
; T = Temperature (Kelvin)
; np = number density of protons (cm^-3)
; v = velocity in reference frame of bulk speed 

common constants_cgs,mp,kb,kev_erg

f1=np*(mp/(2.*!pi*kb*T))^1.5
f2=exp(-1.*mp*(v^2.)/(2.*kb*T))
fm=f1*f2

if fm ne fm then stop
return, fm

end
;------------------------------
function zank_fn,v_plasma,np,del,el,em,ea

common constants_cgs, mp,kb,kev_erg

q=5.
va=sqrt(2*ea*kev_erg/mp)
vl=sqrt(2*el*kev_erg/mp)
vm=sqrt(2*em*kev_erg/mp)
c=np*(del-3.)*(q-3.)/(4*!pi*vl^3.*(q-del))
if v_plasma lt vm then fz=c*((v_plasma/vl)^(-del)-(v_plasma/vl)^(-q))
if v_plasma ge vm then fz=c*((vm/vl)^(q-del)-1.)*(v_plasma/vl)^(-q)

if vl eq 0. then fz=0.
if fz lt 0 then fz=0.

return, fz
end
