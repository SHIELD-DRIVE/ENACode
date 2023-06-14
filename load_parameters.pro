pro load_parameters

;program to define any variables into common blocks to later be called

;grid and data arrays are defined in load_MODELER programs
;all other common blocks are defined here

common grid_params,dtheta,dphi,dr,nphi,ntheta,phi_first,phi_last,theta_first,theta_last,ri_opher,rf_opher,ri_heerikhuisen,rf_heerikhuisen,r_lower,r_upper,r_inner,moscow_regs
dtheta = 6.;4. ; latitude resolution in degrees 
dphi = 6.;4. ;longitude resolution in degrees 
dr = 2. ;radial resolution in AU 
phi_first=0. + dphi/2.;1st phi
phi_last=360.-dphi/2. ;360.- dphi/2.;360.-dphi; last phi
nphi=1.+(phi_last-phi_first)/dphi; number of phi bins
theta_first=0. + dtheta/2.;1st theta
theta_last=180.-dtheta/2.;180. - dtheta/2.;180. - dtheta;last theta
ntheta=1.+(theta_last-theta_first)/dtheta; number of theta bins
ri_opher = 30. ;8. ;AU, smallest radius in grid
rf_opher = 1500. ;1496. ;800.; AU, ;largest radius in grid
ri_heerikhuisen=9.; 
rf_heerikhuisen=1197.;
r_lower=2. ; start point for integration in AU, 30 is inner boundary - MZK
r_upper=800. ; end point for integration (i.e. 252, 402, 804, 1500) - MZK
r_inner=2. ; inner boundary of solution being used
moscow_regs=0 ; 1 = if using region numbers from Moscow then convert to our region numbers
;-----------------------------------------------------------------
common kappa_value,kappa,kappaPUI,kappaRef,kappaISM,kappaOH
kappa= 1.00; kappa value for Kappa distribution (if used)
kappaPUI= 1.00;1.60;10.00; 10.00;1.00; 4.65;1.0; Kappa value for transmitted PUIs - MZK
kappaRef= 1.00;1.60;1.57; 1.57;1.00; Kappa value for reflected PUIs
kappaISM=1.00; kappa value for ISM - MZK
kappaOH= 1.60;kappa value for open heliosheath - MZK
;----------------------------------------------------------------
common energybands,num_ebands,arr_eband_indexes
num_ebands =1; 14.;1. ;3.;14. ; number of energy bands to gp up to, maximum is 14
arr_eband_indexes = [13,9,10,11,12,13,14,15,16,17,20];9,10,11,12,13];0-7 (IBEX-Lo), 8-13 (IBEX-Hi), 14-17 (INCA), 18-20 (LECP);which energy band indexes to run for 
;----------------------------------------------------------------
common radial_limits, radial_range
radial_range = 2 ;0 for TS to HP, 1 for HP to 800 AU, and 2 for TS to 800
;----------------------------------------------------------------
common distributions,dist_names,num_dist
dist_names =['Kappa'+string(kappa,format='(f4.2)')];['Maxwellian'];['shell'];['Kappa'+string(kappa,format='(f4.2)')];
num_dist= 1. ; Default to Kappa distribution since in code have set that kappa < 1.5 means Maxwellian
;----------------------------------------------------------------
common model_names,models,num_models,pui_model,vasyliunas,ribbon,gdf,hinterp,ntr_frac,nref_frac,ttr_frac,tref_frac,colde,esw_frac,etr_frac,eref_frac,moscow,moscow_temp,giacalone
models =['Opher_multiion_withneutrals'];['Opher_multiion_withneutrals'];['Opher_withneutrals'];['Opher'];;; model used from code
num_models = 1.
pui_model=['Malama'];['Malama'];['Zirnstein'];['Other'] ; PUI model (Malama et al. 2006; Zirnsteint et al. 2017)
vasyliunas=0; 0; 0 to not include vasyliunas distribution; PUIs varying from nose to tail
giacalone=1; to use empirical exponential functions for ratios instead of above
ribbon=0; 0; 0 to not include ribbon
gdf=0; 0 for model results, 1 for IBEX observations, 2 for INCA observations 
moscow=0; reading in Moscow ENA model data for plotting -- ENA results from Moscow
moscow_temp=0; 1=if using Moscow MHD model, do not divide by 2 for temp
hinterp=0 ; 0 to not interpolate neutrals, 1 to interpolate neutrals in artificially enhanced locations
ntr_frac=0.23;0.151;0.23;0.151;0.25;0.151;0.168;0.151;0.208;0.150;0.208;0.150;0.208; for Malama method
nref_frac=0.03;0.013;0.03;0.013;0.03;0.013;0.06;0.013;0.04;0.015;0.018;0.013;0.04; for Malama method
ttr_frac=2.17;3.312;2.17;3.312;1.;3.312;2.;3.312;2.21;2.6;3.34;2.21; for Malama method
tref_frac=14.33;35.273;14.33;35.273;21.67;35.273;10.575;35.273;10.575;25.56*1.38;;27.78;35.37;10.575; for Malama method
colde=1; 1 to use cold electrons
esw_frac=0.04 ; for Zirnstein method
etr_frac=0.50 ; for Zirnstein method
eref_frac=0.46 ; for Zirnstein method
;-----------------------------------------------------------------
common include_survival, flag_survival, extinction
flag_survival = 1; 1 ; 0 to not inlude survival prob, for making flux maps at outer boundary
extinction = 1; 1 ; 0 to not include extinction
;-----------------------------------------------------------------
common reading_data,Nlines
Nlines = 3375000;14399 ;7199 ;705599 ;61244252; Nodes-1
;-----------------------------------------------------------------
common path_names,countrate_path,averagefluxplots_path,countrateplots_path,count_input_path,count_output_path
countrate_path='secondary_files/'
averagefluxplots_path='secondary_files/'
count_input_path='input/'
count_output_path='output/'
;-----------------------------------------------------------------
common data_path_names,original_path, opher_version,opher_data_start,heerikhuisen_nh_version,heerikhuisen_version,opher_withneutrals_version, opher_twoplasmas_version, boundaries_version
original_path = 'secondary_files/'; where raw opher data kept
opher_version = '3d_interpolated/3d_Interpolated_SI_400K.dat'
opher_data_start =16 ;34; 16; 16 ;36;   
opher_withneutrals_version = '3d_interpolated/3d_Interpolated_SI_400K.dat';'gridtest1.dat';'Monopole.dat';'tmp1.dat'
opher_twoplasmas_version = 'PLACEHOLDER.dat'
heerikhuisen_nh_version = 'nH_a45_eclip_nathan.ascii'
heerikhuisen_version = 'a45_eclip_plasma.dat'
boundaries_version = '3d_interpolated/3d_Interpolated_SI_400K_withB.dat';'Monopole_withB.dat'
;------------------------------------------------------------------
common secondary_path_names,secondary_path, opher_secondary_file, heerikhuisen_secondary_file, heerikhuisen_secondarynh_file, opher_withneutrals_secondary_file, opher_withneutrals_secondarynh_file, opher_twoplasmas_secondaryp1_file,opher_twoplasmas_secondaryp2_file, opher_twoplasmas_secondarynh_file, opher_withneutrals_secondarystream_file,opher_withneutrals_secondary_file1, opher_withneutrals_secondarynh_file1,opher_withneutrals_secondary_file2, opher_withneutrals_secondarynh_file2,opher_withneutrals_secondaryinterp_file,input_path,output_path
secondary_path = 'secondary_files/'
input_path='input/'
output_path='output/'
opher_secondary_file = 'opherplasmafinish.dat'
opher_withneutrals_secondary_file = 'opherwithneutralsplasmafinish.dat' ; plasma file
opher_withneutrals_secondary_file1 = 'opherwithneutralsplasmafinish1.dat'
opher_withneutrals_secondary_file2 = 'opherwithneutralsplasmafinish2.dat'
opher_withneutrals_secondarynh_file ='opherwithneutralsnhfinish.dat' ; neutral file
opher_withneutrals_secondarynh_file1 ='opherwithneutralsnhfinish1.dat'
opher_withneutrals_secondarynh_file2 ='opherwithneutralsnhfinish2.dat'
opher_withneutrals_secondaryinterp_file='opherinterpolatednhinput.dat'
opher_withneutrals_secondarystream_file='Streamlines/streamlines_kmhd_BISMtest_Alpha10_local_IB2AU_2AU6deg.dat';'Streamlines/streamlines_moscow_nohe_2AU.dat';'Streamlines/streamlines_kmhd_moscow_n1200_IB2AU_2AU.dat';'Streamlines/streamlines_moscow_IA2020Cart2AU_IB2AU_2AU_2.dat';'Streamlines/streamlines_kmhd_2014_hires_2AU.dat';'Streamlines/streamlines_kmhd_LoB_uniform_n50900_2AU_smoothed2.dat';'Streamlines/streamlines_kmhd_nobsw_hires_2AU.dat';'Streamlines/streamlines_kmhd_withbsw_hires_2AU.dat'; for extinction
heerikhuisen_secondary_file = 'heerikhuisenplasmafinishcb.dat'
heerikhuisen_secondarynh_file = 'heerikhuisennhfinish.dat'
;-------------------------------------------------------------------
common plot_params,latmin,latmax,lonmin,lonmax,map_dlat,map_dlon,map_numlatlabels,map_lons,map_lonnames,onerange_flux_min, onerange_flux_max, onerange_cr_min, onerange_cr_max
latmin = -90.;in degrees, in IDL space with nose at 0
latmax = 90.
lonmin = -180.
lonmax = 180.
map_dlat=30.
map_dlon=45.
map_numlatlabels = 8
map_lons =  [-180,-135,-90,-45,0,45,90,135,180];IDL assumes low to high
map_lonnames = [75,30,345,300,255,210,165,120,75]; ;from high to low for sky map
onerange_flux_min = [10000.,10000.,10000.,10000.,10000.] ;[0.,0.,0.,0.,0.];ARRAY of min for each energyband for manual forcing of min and max of plot, good to use for special cases, such as plotting small portion of sky
onerange_flux_max = [400000.,400000.,400000.,400000.,400000.];[130.,130.,130.,130.,130.] ;ARRAY of max for each energyband, should be same size as num_ebands
onerange_cr_min = [0., 0., 0., 0., 0.]; Size of this array should be the same as the num_ebands
onerange_cr_max = [2., 2., 2., 2., 2.]; Size of this array should be the same as the num_ebands
;-----------------------------------------------------------------
common plot_flags, onerange_flag, manual_onerange_flag, plot_log, plot_limited, manual_lim, manual_maxset, manual_minset, eps_plot, plot_scl,tail, spec_min, spec_max, pres_min, pres_max, write_ecl,legend_style
onerange_flag = 0; if 1 then find max and min for each energy band over all modelers and distributions for comparison, set to 0 if plotting small section of sky and either use manual range or allow pro to calc. max of truncated plot
manual_onerange_flag =0.;  set to 1 if you want to specify min and max (in plot_params), best to use when only using 1 or 2 energy channels
plot_log = 0 ;0. ; if 1 then plot intensity on log scale
plot_limited = 0. ;set to one for now
manual_lim=1. ; 1 for manual color bar limits, 0 for default color bar limits - MZK
manual_maxset=[12.,265.,140.,60.,30.,12.,0.90,0.08,0.025,0.009,5e-4];[800.,400.,140.,60.,18.][265.,140.,60.,30.,12.][120.,60.,25.,10.,6.]
manual_minset=[0.1,0,0,0,0,0,0,0,0,0] ;0.195, min color bar value when manual_lim=1 - MZK
eps_plot=1; if 1, then plot in eps format, otherwise plot in jpg
plot_scl=1.0;1.8 ; scale flux in ENA maps
tail=0; if 1, center map on tail; otherwise, center map on nose
spec_min=[0.]; min color bar value for spectral slope plot
spec_max=[6.]; max color bar value for spectral slope plot
pres_min=[0.]; min color bar value for pressure plot
pres_max=[30.]; max color bar value for pressure plot
write_ecl=0 ; via plotter_eps.pro - write data file in ecliptic coords
legend_style=1 ; number format for color bar legend: 0=integer, 1=float, 2=exponent (small numbers)
;-----------------------------------------------------------------
common log_def, log_num
log_num = 10.
;-------------------------------------------------------------------
common constants_cgs, mp,kb,kev_erg,AU_cm
mp=1.66D*10.0^(-24D);mass of a proton in grams
kb=1.3806503D*10.0^(-16D)  ; cm^2 g s^-2 K^-1
kev_erg = 1.60217646D*10^(-9D) ;conversion from keV to ergs
AU_cm = 1.49598D*10^(13D)
;------------------------------------------------------------------
common sensor_params, hi,lo, wt9, wt10, wt11, wt12, wt13, ultra
;Loads information about IBEX HI and Lo into a commonblock (min energies,
;max energies, central energies, and geometric factors)
;From Sources Document 

Energyultra={ultraenergies, min:fltarr(8), max:fltarr(8), central:fltarr(8)}
Energyhi =  {hienergies, min:fltarr(6), max:fltarr(6), central:fltarr(6)}
Energylo =  {loenergies, min:fltarr(8), max:fltarr(8), central:fltarr(8)}
Ultra= {Ultra_Sensor, E:Energyultra, G:fltarr(8), Noise:fltarr(8) }
Hi = {Hi_Sensor, E:Energyhi, G:fltarr(6), Noise:fltarr(6) }
lo = {lo_Sensor, E:Energylo, G:fltarr(8), Noise:fltarr(8) }

Ultra.E.min     = [5.2, 13.5, 24., 35., 28., 60.,70.];43.0, 80.0, 137.0] ; INCA + V2 energies
Ultra.E.max     = [13.5, 24., 35., 55., 43., 80.,90.];80.0,137.0, 215.0]
Ultra.E.central = [8.37854, 18.0, 28.9828, 43.8748, 33.4664, 70.,80.];56.5685, 104.690,171.625]
Hi.E.min     = [0.38, 0.52, 0.84, 1.36, 1.99, 3.13] ; HW,triple coincidence energy values from Table  6 Funsten 2009
Hi.E.max    =  [0.59, 0.95, 1.55, 2.5,  3.75, 6.0]
Hi.E.central = [0.45, 0.71, 1.11, 1.74, 2.73, 4.29]
Lo.E.min     = [0.009, 0.016,0.031,0.061,0.118, 0.271,0.545,1.142]
Lo.E.max     = [0.020, 0.038,0.073,0.143,0.276, 0.631,1.27, 2.66 ]
Lo.E.central =  [0.014,0.027,0.052,0.102,0.197, 0.451,0.908,1.903] ;from Fuselier et al., 2009

Ultra.G = [5.0E-3, 5.5E-3, 6.0E-3, 6.5E-3, 7.0E-3, 7.5E-3, 8.0E-3, 8.5E-3]
Hi.G = [1.3E-4, 4.1E-4, 7.5E-4,1.3E-3, 2.4E-3, 4.5E-3]; funsten et al 2009, triples only - fixed by MZK
Lo.G = [2.6E-5, 5.3E-5, 8.1E-5, 9.1E-5, 9.0E-5, 1.0E-4, 1.9E-4, 2.7E-4] ; from Fuselier et al., 2009, triples only

;Noise rates for each energy channel, s^-1
Ultra.Noise = [0., 0., 0., 0., 0., 0., 0., 0.]
Hi.Noise = [ 3.8e-2, 6.8e-5, 1.2e-4, 1.9e-4, 2.6e-4, 4.2e-4]
Lo.Noise = [ 4.7e-4, 8.9e-4, 1.7e-3, 3.2e-3, 6.3e-3, 1.9e-2, 0. , 0. ]

Wt9=[0.0013277909,0.0016871014,0.0019651474,0.0022037048,0.0022894956,0.0022060268,0.0020821493,0.001934731,0.0016461171,0.0012819526]
Wt10=[0.0008472022,0.0011271987,0.0013325584,0.0014113985,0.0013964271,0.0013630371,0.0012635717,0.0011411818,0.0009836919,0.000787456]
Wt11=[0.0005482951,0.0007248594,0.0008530854,0.0008994173,0.0008982827,0.000832246,0.0007757703,0.000709315,0.0006125453,0.0005138144]
Wt12=[0.0003588889,0.0004682034,0.0005517594,0.0005910363,0.0005914153,0.0005547477,0.000513061,0.0004612037,0.0004018451,0.0003388693]
Wt13=[0.0002324185,0.0002991092,0.000347635,0.0003659699,0.0003612539,0.000345767,0.0003205023,0.0002890332,0.0002574077,0.000213177]
;--------------------------------------------------------------------

return
end
