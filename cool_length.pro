pro cool_length, cool

  ;num=1000
  ;energy=findgen(num)/10.+0.1

  ;cool=cool*0.1/0.155
  
  num=12
  energy=[0.1,0.7,1.1,1.7,2.7,10,20,30,40,50,75,100]
  
  vena=energy_velocity(energy)

  a1=4.15
  a2=0.531
  a3=67.3

  sig1=(a1-a2*alog(energy))^2
  sig2=(1-exp(-a3/energy))^(4.5)
  sigma=sig1*sig2
  sigma=sigma*10D^(-16D)

  ext=make_array(num,736)
  for i=0,num-1 do ext(i,*)=cool[12,3,*]*vena[i]*sigma[i]

  lc=make_array(num)

  for j=0,num-1 do begin
  depth=where(abs(1/exp(1)+ext[j,*]) eq min(abs(1/exp(1)+ext[j,*])),dct)
  lc[j]=min(depth)
  endfor                                                     

  p=plot(energy,lc*2., /xlog,/ylog, xtitle="Energy [keV]", ytitle="Cooling Length [AU]", title="Variation of Cooling Length with Energy (nH=0.18 cc)", yrange=[1,1e3],xrange=[0.1,500])
  ;p=plot(energy,sigma, /xlog, /ylog, xtitle="Energy [keV]", ytitle="Cross-section [cm^2]", title="Variation of Cross-section with Energy",yrange=[1e-18,1e-14])
  ;p.save, "cooling_length.png"
  ;p.save, "cross_section.png"

  stop
 end
