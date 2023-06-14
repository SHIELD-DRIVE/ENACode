pro load_opher
;Creates common block of r, theta, and phi grid indexes based on size
;of Merav's grid

common grid_opher,gridopher
common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen
common data_opher,dataopher,sdataopher
common nhdata_opher,nhopher,snhopher ;used when neutrals

nr=1.+floor((rf_opher-ri_opher)/dr); number of r bins, added floor

separations = {indivsep, first: 1., last:10., delta:1., num:10.}
gridopher= {allsep, phi:separations, theta:separations, r:separations}

gridopher.phi.first=phi_first
gridopher.phi.last=phi_last
gridopher.phi.delta=dphi
gridopher.phi.num=nphi
gridopher.theta.first=theta_first
gridopher.theta.last=theta_last
gridopher.theta.delta=dtheta
gridopher.theta.num=ntheta
gridopher.r.first=ri_opher
gridopher.r.last=rf_opher
gridopher.r.delta=dr
gridopher.r.num=nr

;creates structure for processing secondary data
other=dblarr(ntheta,nphi,nr)
den=dblarr(3,ntheta,nphi,nr)
neutral=dblarr(ntheta,nphi,nr,4)
dataopher={plasmaopher, density:den, temp:other,vr:other,vphi:other,vtheta:other,vt:other,theta:other,phi:other,r:other}
nhopher = {neutralopher, density:other}
;nhopher = {neutralopher,density:neutral, vr:neutral, vt:neutral}
;creates structure for loading secondary
sdataopher={splasmaopher, density:other, temp:other,vr:other,vphi:other,vtheta:other,vt:other,theta:other,phi:other,r:other}
snhopher = {sneutralopher,density:other, temp:other}
;snhopher = {sneutralopher,density:neutral, vr:other, vt:other, temp:other} ; for neutral pops -MZK

end

