#!MC 1410
# Step 3 of the single-ion input preparation with the streamlines traced by ena_trace
# (configs/prep_single_ion.toml).
# Starts from: the spherical grid with the interpolated neutrals and plasma (after step 2).
# Writes, in this order: |NEUTRALFILE| (8 columns), |PLASMAFILE| (12), |BOUNDARYFILE| (10, with the
#       region definitions), and |STREAMFIELD|: the spherical grid with the 19 streamline variables,
#       which ena_trace traces the streamlines through (instead of Tecplot's streamtraces, which are
#       not reproducible).
# The same as 3_calculate_files.mcr up to the streamtraces; the 1800 $!Streamtrace Add commands and
#       the extraction are replaced by the write of the field. The seeds are in
#       configs/streamline_seeds_2AU.txt, read by ena_trace.
$!AlterData 
  Equation = '{r}=sqrt({X AU}**2+{Y AU}**2+{Z AU}**2)'
$!AlterData 
  Equation = '{theta}=acos({Z AU}/({r}+1e-4))*57.2958'
$!AlterData 
  Equation = '{phi}=atan2({Y AU},{X AU})*57.2958'
$!AlterData 
  Equation = '{phi}=if({phi} < 0, 360 - abs({phi}), {phi})'
$!AlterData 
  Equation = '{v1}={nH}'
$!AlterData 
  Equation = '{v2}={TH}'
$!AlterData 
  Equation = '{v3}=k-1'
$!AlterData 
  Equation = '{v4}=i-1'
$!AlterData 
  Equation = '{v5}=j-1'
$!AlterData 
  Equation = '{v6}={theta}'
$!AlterData 
  Equation = '{v7}={phi}'
$!AlterData 
  Equation = '{v8}={r}'
$!WriteDataSet  "|NEUTRALFILE|"
  IncludeText = No
  IncludeGeom = No
  IncludeCustomLabels = No
  IncludeDataShareLinkage = Yes
  VarPositionList =  [44-51]
  Binary = No
  UsePointFormat = Yes
  Precision = 9
  TecplotVersionToWrite = TecplotCurrent
$!DeleteVars  [41-51]
$!AlterData 
  Equation = '{r}=sqrt({X AU}**2+{Y AU}**2+{Z AU}**2)'
$!AlterData 
  Equation = '{theta}=acos({Z AU}/({r}+1e-4))*57.2958'
$!AlterData 
  Equation = '{phi}=atan2({Y AU},{X AU})*57.2958'
$!AlterData 
  Equation = '{phi}=if({phi} < 0, 360 - abs({phi}), {phi})'
$!AlterData 
  Equation = '{v1}=v4'
$!AlterData 
  Equation = '{v2}=v11/(v4*1.38e-16)'
$!AlterData 
  Equation = '{v3}=({U_x km/s}*{X AU}+{U_y km/s}*{Y AU}+{U_z km/s}*{Z AU})/({r}+1e-4)'
$!AlterData 
  Equation = '{v4}={r}*({U_x km/s}*{Y AU}-{U_y km/s}*{X AU})/({X AU}**2+{Y AU}**2+1e-4)'
$!AlterData 
  Equation = '{v5}={r}*sin({theta}*(0.0174533))*({Z AU}*({X AU}*{U_x km/s}+{Y AU}*{U_y km/s})-({X AU}**2+{Y AU}**2)*{U_Z km/s})/(({X AU}**2+{Y AU}**2+{Z AU}**2)*sqrt({X AU}**2+{Y AU}**2)+1e-4)'
$!AlterData 
  Equation = '{v6}=sqrt(({U_x km/s}**2+{U_y km/s}**2+{U_z km/s}**2)-{v3}**2)'
$!AlterData 
  Equation = '{v7}=k-1'
$!AlterData 
  Equation = '{v8}=i-1'
$!AlterData 
  Equation = '{v9}=j-1'
$!AlterData 
  Equation = '{v10}={theta}'
$!AlterData 
  Equation = '{v11}={phi}'
$!AlterData 
  Equation = '{v12}={r}'
$!WriteDataSet  "|PLASMAFILE|"
  IncludeText = No
  IncludeGeom = No
  IncludeCustomLabels = No
  IncludeDataShareLinkage = Yes
  VarPositionList =  [44-55]
  Binary = No
  UsePointFormat = Yes
  Precision = 9
  TecplotVersionToWrite = TecplotCurrent
$!DeleteVars  [41-55]
$!AlterData 
  Equation = '{X AU}=v1'
$!AlterData 
  Equation = '{Y AU}=v2'
$!AlterData 
  Equation = '{Z AU}=v3'
$!AlterData 
  Equation = '{rho}=v4*1.67e-24'
$!AlterData 
  Equation = '{Cs^2}=5*v11/(({rho}+1e-30)*3)'
$!AlterData 
  Equation = '{U^2}=v5**2+v6**2+v7**2'
$!AlterData 
  Equation = '{U^2}={U^2}*1e10'
$!AlterData 
  Equation = '{M^2}={U^2}/({Cs^2}+0.0001)'
$!AlterData 
  Equation = '{T [K]}=v11/(1.38e-16*(v4+1e-30))'
$!AlterData 
  Equation = '{lnT}=log({T [K]}+0.0001)'
$!AlterData 
  Equation = '{r}=({X AU}**2+{Y AU}**2+{Z AU}**2)**(0.5)'
$!AlterData 
  Equation = '{region}=if({M^2} > 1, if (v40 > 0, 4, 5), 5)'
$!AlterData 
  Equation = '{region}=if({M^2} < 1,3,{region})'
$!AlterData 
  Equation = '{region}=if(v40<0, if({U^2}<1e14, 2, {region}), {region})'
$!AlterData 
  Equation = '{region}=if({M^2} > 1, if({U^2} < 1e14, if({r} > 400, 1, {region}), {region}), {region})'
$!AlterData 
  Equation = '{theta}=acos({Z AU}/({r}+1e-10))*57.29577'
$!AlterData 
  Equation = '{phi}=atan2({Y AU},{X AU})*57.29577'
$!AlterData 
  Equation = '{phi}=if({phi}<0, 360-abs({phi}), {phi})'
$!AlterData 
  Equation = '{v1}={region}'
$!AlterData 
  Equation = '{v2}=k-1'
$!AlterData 
  Equation = '{v3}=i-1'
$!AlterData 
  Equation = '{v4}=j-1'
$!AlterData 
  Equation = '{v5}={theta}'
$!AlterData 
  Equation = '{v6}={phi}'
$!AlterData 
  Equation = '{v7}={r}'
$!AlterData 
  Equation = '{v8}=v8'
$!AlterData 
  Equation = '{v9}=v9'
$!AlterData 
  Equation = '{v10}=v10'
$!WriteDataSet  "|BOUNDARYFILE|"
  IncludeText = No
  IncludeGeom = No
  IncludeCustomLabels = No
  IncludeDataShareLinkage = Yes
  VarPositionList =  [51-60]
  Binary = No
  UsePointFormat = Yes
  Precision = 9
  TecplotVersionToWrite = TecplotCurrent
$!DeleteVars  [41-60]
$!AlterData 
  Equation = '{r}=sqrt({X AU}**2+{Y AU}**2+{Z AU}**2)'
$!AlterData 
  Equation = '{theta}=acos({Z AU}/({r}+1e-4))*57.2958'
$!AlterData 
  Equation = '{phi}=atan2({Y AU},{X AU})*57.2958'
$!AlterData 
  Equation = '{phi}=if({phi} < 0, 360 - abs({phi}), {phi})'
$!AlterData 
  Equation = '{np}=v4'
$!AlterData 
  Equation = '{tp}=v11/(v4*1.38e-16)'
$!AlterData 
  Equation = '{vr}=({U_x km/s}*{X AU}+{U_y km/s}*{Y AU}+{U_z km/s}*{Z AU})/({r}+1e-4)'
$!AlterData 
  Equation = '{vtheta}={r}*({U_x km/s}*{Y AU}-{U_y km/s}*{X AU})/({X AU}**2+{Y AU}**2+1e-4)'
$!AlterData 
  Equation = '{vphi}={r}*sin({theta}*(0.0174533))*({Z AU}*({X AU}*{U_x km/s}+{Y AU}*{U_y km/s})-({X AU}**2+{Y AU}**2)*{U_Z km/s})/(({X AU}**2+{Y AU}**2+{Z AU}**2)*sqrt({X AU}**2+{Y AU}**2)+1e-4)'
$!AlterData 
  Equation = '{vt}=sqrt(({U_x km/s}**2+{U_y km/s}**2+{U_z km/s}**2)-{vr}**2)'
$!AlterData 
  Equation = '{theta0}={theta}'
$!AlterData 
  Equation = '{phi0}={phi}'
$!AlterData 
  Equation = '{r0}={r}'
# The grid (zone 1) with the variables the streamline file has (as Tecplot's streamtraces would
# carry them): X, Y, Z, U_x, U_y, U_z, nH, UxH, UyH, UzH, np, tp, vr, vtheta, vphi, vt, theta0,
# phi0, r0, in point format (ena_trace reads it).
$!WriteDataSet  "|STREAMFIELD|"
  IncludeText = No
  IncludeGeom = No
  IncludeCustomLabels = No
  IncludeDataShareLinkage = Yes
  ZoneList =  [1]
  VarPositionList =  [1-3,5-7,35-38,44-52]
  Binary = No
  UsePointFormat = Yes
  Precision = 9
  TecplotVersionToWrite = TecplotCurrent
