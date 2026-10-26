#!MC 1410
# Step 3 of the single-ion input preparation (tools/tecplot_prep.py).
# Starts from: the spherical grid with the interpolated neutrals and plasma (after step 2).
# Writes, in this order: |NEUTRALFILE| (8 columns), |PLASMAFILE| (12), |BOUNDARYFILE| (10, with the
#       region definitions), and |STREAMFILE| (the 1800 streamtraces seeded at 2 AU, extracted;
#       19 variables).
# Generalised from calculating_files.mcr: the four output file names were changed, and the zone list
#       of the streamline file is [2-|NUMZONES|] (all extracted streamtraces, however many a solution
#       gives) instead of the fixed [2-1583] of the recorded solution. The first, accidental copy of
#       the region calculation was removed; the plasma-file variables 41-55 are
#       deleted before the regions are computed, so the variable numbers written are unchanged.
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
$!StreamtraceLayers Show = Yes
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.209056928753853
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.415823370218277
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.6180340051651
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.813473284244537
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.17557048797607
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33826124668121
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.48628962039948
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6180340051651
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73205077648163
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.82709085941315
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.90211308002472
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.95629525184631
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.98904383182526
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 2
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.98904383182526
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.95629525184631
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.90211308002472
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.82709085941315
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73205077648163
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6180340051651
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.48628962039948
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33826124668121
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.17557048797607
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.813473284244537
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.6180340051651
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.415823370218277
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.209056928753853
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 2.44929370517034E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.209056928753853
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.415823370218277
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.6180340051651
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.813473284244537
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.17557048797607
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33826124668121
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.48628962039948
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6180340051651
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73205077648163
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.82709085941315
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.90211308002472
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.95629525184631
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.98904383182526
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -2
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.98904383182526
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.95629525184631
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.90211308002472
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.82709085941315
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73205077648163
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6180340051651
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.48628962039948
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33826124668121
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.17557048797607
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.813473284244537
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.6180340051651
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.415823370218277
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.209056928753853
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = 0.0218524001538754
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = 0.0434653796255589
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.614648342132568
    Y = 0.0646021440625191
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = 0.0850311145186424
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = 0.104528464376926
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.16913056373596
    Y = 0.122880578041077
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33093011379242
    Y = 0.139886394143105
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = 0.1553595662117
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.60917019844055
    Y = 0.169130608439445
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.72256243228912
    Y = 0.181048601865768
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.81708192825317
    Y = 0.190983012318611
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.89169311523438
    Y = 0.19882495701313
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.94557845592499
    Y = 0.204488530755043
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.97814762592316
    Y = 0.207911685109138
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.98904383182526
    Y = 0.209056928753853
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.97814762592316
    Y = 0.207911685109138
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.94557845592499
    Y = 0.204488530755043
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.89169311523438
    Y = 0.19882495701313
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.81708192825317
    Y = 0.190983012318611
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.72256243228912
    Y = 0.181048601865768
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.60917019844055
    Y = 0.169130608439445
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = 0.1553595662117
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33093011379242
    Y = 0.139886394143105
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.16913056373596
    Y = 0.122880578041077
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = 0.104528464376926
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = 0.0850311145186424
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.614648342132568
    Y = 0.0646021440625191
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = 0.0434653796255589
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = 0.0218524001538754
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 2.43587617389024E-16
    Y = 2.56020903001859E-17
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = -0.0218524001538754
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = -0.0434653796255589
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.614648342132568
    Y = -0.0646021440625191
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = -0.0850311145186424
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = -0.104528464376926
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.16913056373596
    Y = -0.122880578041077
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33093011379242
    Y = -0.139886394143105
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = -0.1553595662117
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.60917019844055
    Y = -0.169130608439445
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.72256243228912
    Y = -0.181048601865768
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.81708192825317
    Y = -0.190983012318611
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.89169311523438
    Y = -0.19882495701313
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.94557845592499
    Y = -0.204488530755043
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.97814762592316
    Y = -0.207911685109138
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.98904383182526
    Y = -0.209056928753853
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.97814762592316
    Y = -0.207911685109138
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.94557845592499
    Y = -0.204488530755043
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.89169311523438
    Y = -0.19882495701313
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.81708192825317
    Y = -0.190983012318611
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.72256243228912
    Y = -0.181048601865768
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.60917019844055
    Y = -0.169130608439445
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = -0.1553595662117
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33093011379242
    Y = -0.139886394143105
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.16913056373596
    Y = -0.122880578041077
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = -0.104528464376926
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = -0.0850311145186424
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.614648342132568
    Y = -0.0646021440625191
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = -0.0434653796255589
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = -0.0218524001538754
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.204488530755043
    Y = 0.0434653796255589
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = 0.0864545404911041
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = 0.128496497869492
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.795696914196014
    Y = 0.169130608439445
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.978147625923157
    Y = 0.207911685109138
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.14988148212433
    Y = 0.24441485106945
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = 0.278240144252777
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.4538106918335
    Y = 0.30901700258255
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58267605304718
    Y = 0.336408168077469
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6942013502121
    Y = 0.360113620758057
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.78716456890106
    Y = 0.37987357378006
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.86054730415344
    Y = 0.395471543073654
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.91354548931122
    Y = 0.406736642122269
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.94557845592499
    Y = 0.413545459508896
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.95629525184631
    Y = 0.415823370218277
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.94557845592499
    Y = 0.413545459508896
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.91354548931122
    Y = 0.406736642122269
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.86054730415344
    Y = 0.395471543073654
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.78716456890106
    Y = 0.37987357378006
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6942013502121
    Y = 0.360113620758057
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58267605304718
    Y = 0.336408168077469
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.4538106918335
    Y = 0.30901700258255
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = 0.278240144252777
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.14988148212433
    Y = 0.24441485106945
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.978147625923157
    Y = 0.207911685109138
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.795696914196014
    Y = 0.169130608439445
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = 0.128496497869492
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = 0.0864545404911041
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.204488530755043
    Y = 0.0434653796255589
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 2.39577075202452E-16
    Y = 5.09236770372102E-17
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.204488530755043
    Y = -0.0434653796255589
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = -0.0864545404911041
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = -0.128496497869492
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.795696914196014
    Y = -0.169130608439445
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.978147625923157
    Y = -0.207911685109138
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.14988148212433
    Y = -0.24441485106945
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = -0.278240144252777
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.4538106918335
    Y = -0.30901700258255
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58267605304718
    Y = -0.336408168077469
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6942013502121
    Y = -0.360113620758057
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.78716456890106
    Y = -0.37987357378006
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.86054730415344
    Y = -0.395471543073654
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.91354548931122
    Y = -0.406736642122269
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.94557845592499
    Y = -0.413545459508896
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.95629525184631
    Y = -0.415823370218277
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.94557845592499
    Y = -0.413545459508896
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.91354548931122
    Y = -0.406736642122269
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.86054730415344
    Y = -0.395471543073654
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.78716456890106
    Y = -0.37987357378006
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6942013502121
    Y = -0.360113620758057
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58267605304718
    Y = -0.336408168077469
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.4538106918335
    Y = -0.30901700258255
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = -0.278240144252777
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.14988148212433
    Y = -0.24441485106945
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.978147625923157
    Y = -0.207911685109138
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.795696914196014
    Y = -0.169130608439445
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = -0.128496497869492
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = -0.0864545404911041
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.204488530755043
    Y = -0.0434653796255589
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.19882495701313
    Y = 0.0646021440625191
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.395471543073654
    Y = 0.128496497869492
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = 0.190983012318611
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.77365905046463
    Y = 0.251377075910568
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = 0.30901700258255
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.1180340051651
    Y = 0.363271266222
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.27276206016541
    Y = 0.413545459508896
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.41354548931122
    Y = 0.459288746118546
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.53884172439575
    Y = 0.5
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.64727818965912
    Y = 0.53523313999176
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73766672611237
    Y = 0.564602136611938
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.80901694297791
    Y = 0.587785243988037
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.86054730415344
    Y = 0.604528486728668
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.89169311523438
    Y = 0.614648342132568
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.90211308002472
    Y = 0.6180340051651
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.89169311523438
    Y = 0.614648342132568
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.86054730415344
    Y = 0.604528486728668
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.80901694297791
    Y = 0.587785243988037
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73766672611237
    Y = 0.564602136611938
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.64727818965912
    Y = 0.53523313999176
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.53884172439575
    Y = 0.5
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.41354548931122
    Y = 0.459288746118546
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.27276206016541
    Y = 0.413545459508896
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.1180340051651
    Y = 0.363271266222
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = 0.30901700258255
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.77365905046463
    Y = 0.251377075910568
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = 0.190983012318611
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.395471543073654
    Y = 0.128496497869492
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.19882495701313
    Y = 0.0646021440625191
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 2.32941657321678E-16
    Y = 7.56873370685124E-17
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.19882495701313
    Y = -0.0646021440625191
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.395471543073654
    Y = -0.128496497869492
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = -0.190983012318611
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.77365905046463
    Y = -0.251377075910568
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = -0.30901700258255
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.1180340051651
    Y = -0.363271266222
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.27276206016541
    Y = -0.413545459508896
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.41354548931122
    Y = -0.459288746118546
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.53884172439575
    Y = -0.5
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.64727818965912
    Y = -0.53523313999176
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73766672611237
    Y = -0.564602136611938
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.80901694297791
    Y = -0.587785243988037
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.86054730415344
    Y = -0.604528486728668
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.89169311523438
    Y = -0.614648342132568
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.90211308002472
    Y = -0.6180340051651
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.89169311523438
    Y = -0.614648342132568
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.86054730415344
    Y = -0.604528486728668
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.80901694297791
    Y = -0.587785243988037
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73766672611237
    Y = -0.564602136611938
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.64727818965912
    Y = -0.53523313999176
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.53884172439575
    Y = -0.5
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.41354548931122
    Y = -0.459288746118546
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.27276206016541
    Y = -0.413545459508896
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.1180340051651
    Y = -0.363271266222
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = -0.30901700258255
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.77365905046463
    Y = -0.251377075910568
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = -0.190983012318611
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.395471543073654
    Y = -0.128496497869492
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.19882495701313
    Y = -0.0646021440625191
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = 0.0850311145186424
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.37987357378006
    Y = 0.169130608439445
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.564602136611938
    Y = 0.251377075910568
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = 0.330869406461716
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.913545429706573
    Y = 0.406736642122269
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.07393705844879
    Y = 0.478147596120834
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.22256243228912
    Y = 0.544319868087769
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.35779321193695
    Y = 0.604528486728668
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = 0.658113718032837
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58230710029602
    Y = 0.7044885158538819
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.66913056373596
    Y = 0.743144810199738
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73766672611237
    Y = 0.77365905046463
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.78716456890106
    Y = 0.795696914196014
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.81708192825317
    Y = 0.80901700258255
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.82709085941315
    Y = 0.813473284244537
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.81708192825317
    Y = 0.80901700258255
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.78716456890106
    Y = 0.795696914196014
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73766672611237
    Y = 0.77365905046463
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.66913056373596
    Y = 0.743144810199738
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58230710029602
    Y = 0.7044885158538819
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = 0.658113718032837
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.35779321193695
    Y = 0.604528486728668
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.22256243228912
    Y = 0.544319868087769
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.07393705844879
    Y = 0.478147596120834
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.913545429706573
    Y = 0.406736642122269
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = 0.330869406461716
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.564602136611938
    Y = 0.251377075910568
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.37987357378006
    Y = 0.169130608439445
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = 0.0850311145186424
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 2.23754102701048E-16
    Y = 9.96217485357301E-17
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = -0.0850311145186424
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.37987357378006
    Y = -0.169130608439445
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.564602136611938
    Y = -0.251377075910568
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = -0.330869406461716
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.913545429706573
    Y = -0.406736642122269
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.07393705844879
    Y = -0.478147596120834
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.22256243228912
    Y = -0.544319868087769
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.35779321193695
    Y = -0.604528486728668
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = -0.658113718032837
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58230710029602
    Y = -0.7044885158538819
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.66913056373596
    Y = -0.743144810199738
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73766672611237
    Y = -0.77365905046463
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.78716456890106
    Y = -0.795696914196014
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.81708192825317
    Y = -0.80901700258255
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.82709085941315
    Y = -0.813473284244537
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.81708192825317
    Y = -0.80901700258255
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.78716456890106
    Y = -0.795696914196014
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73766672611237
    Y = -0.77365905046463
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.66913056373596
    Y = -0.743144810199738
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58230710029602
    Y = -0.7044885158538819
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = -0.658113718032837
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.35779321193695
    Y = -0.604528486728668
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.22256243228912
    Y = -0.544319868087769
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.07393705844879
    Y = -0.478147596120834
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.913545429706573
    Y = -0.406736642122269
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = -0.330869406461716
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.564602136611938
    Y = -0.251377075910568
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.37987357378006
    Y = -0.169130608439445
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = -0.0850311145186424
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.181048601865768
    Y = 0.104528464376926
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.360113620758057
    Y = 0.207911685109138
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.53523313999176
    Y = 0.30901700258255
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.7044885158538819
    Y = 0.406736642122269
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.866025388240814
    Y = 0.5
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.01807391643524
    Y = 0.587785243988037
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.15896821022034
    Y = 0.669130623340607
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.28716456890106
    Y = 0.743144810199738
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.40125858783722
    Y = 0.80901700258255
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.5
    Y = 0.866025388240814
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58230710029602
    Y = 0.913545429706573
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.64727818965912
    Y = 0.95105654001236
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6942013502121
    Y = 0.978147625923157
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.72256243228912
    Y = 0.994521915912628
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73205077648163
    Y = 1
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.72256243228912
    Y = 0.994521915912628
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6942013502121
    Y = 0.978147625923157
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.64727818965912
    Y = 0.95105654001236
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58230710029602
    Y = 0.913545429706573
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.5
    Y = 0.866025388240814
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.40125858783722
    Y = 0.80901700258255
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.28716456890106
    Y = 0.743144810199738
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.15896821022034
    Y = 0.669130623340607
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.01807391643524
    Y = 0.587785243988037
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.866025388240814
    Y = 0.5
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.7044885158538819
    Y = 0.406736642122269
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.53523313999176
    Y = 0.30901700258255
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.360113620758057
    Y = 0.207911685109138
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.181048601865768
    Y = 0.104528464376926
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 2.12115049442606E-16
    Y = 1.22464685258517E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.181048601865768
    Y = -0.104528464376926
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.360113620758057
    Y = -0.207911685109138
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.53523313999176
    Y = -0.30901700258255
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.7044885158538819
    Y = -0.406736642122269
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.866025388240814
    Y = -0.5
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.01807391643524
    Y = -0.587785243988037
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.15896821022034
    Y = -0.669130623340607
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.28716456890106
    Y = -0.743144810199738
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.40125858783722
    Y = -0.80901700258255
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.5
    Y = -0.866025388240814
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58230710029602
    Y = -0.913545429706573
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.64727818965912
    Y = -0.95105654001236
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6942013502121
    Y = -0.978147625923157
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.72256243228912
    Y = -0.994521915912628
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73205077648163
    Y = -1
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.72256243228912
    Y = -0.994521915912628
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6942013502121
    Y = -0.978147625923157
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.64727818965912
    Y = -0.95105654001236
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58230710029602
    Y = -0.913545429706573
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.5
    Y = -0.866025388240814
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.40125858783722
    Y = -0.80901700258255
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.28716456890106
    Y = -0.743144810199738
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.15896821022034
    Y = -0.669130623340607
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.01807391643524
    Y = -0.587785243988037
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.866025388240814
    Y = -0.5
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.7044885158538819
    Y = -0.406736642122269
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.53523313999176
    Y = -0.30901700258255
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.360113620758057
    Y = -0.207911685109138
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.181048601865768
    Y = -0.104528464376926
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = 0.122880578041077
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.336408168077469
    Y = 0.24441485106945
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = 0.363271266222
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.658113718032837
    Y = 0.478147596120834
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = 0.587785243988037
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = 0.69098299741745
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.08267605304718
    Y = 0.786610186100006
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.20243358612061
    Y = 0.873619139194489
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = 0.95105654001236
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.40125858783722
    Y = 1.01807391643524
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = 1.07393705844879
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.53884172439575
    Y = 1.1180340051651
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58267605304718
    Y = 1.14988148212433
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.60917019844055
    Y = 1.16913056373596
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6180340051651
    Y = 1.17557048797607
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.60917019844055
    Y = 1.16913056373596
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58267605304718
    Y = 1.14988148212433
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.53884172439575
    Y = 1.1180340051651
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = 1.07393705844879
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.40125858783722
    Y = 1.01807391643524
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = 0.95105654001236
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.20243358612061
    Y = 0.873619139194489
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.08267605304718
    Y = 0.786610186100006
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = 0.69098299741745
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = 0.587785243988037
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.658113718032837
    Y = 0.478147596120834
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = 0.363271266222
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.336408168077469
    Y = 0.24441485106945
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = 0.122880578041077
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.98152015709584E-16
    Y = 1.43965862235849E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = -0.122880578041077
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.336408168077469
    Y = -0.24441485106945
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = -0.363271266222
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.658113718032837
    Y = -0.478147596120834
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = -0.587785243988037
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = -0.69098299741745
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.08267605304718
    Y = -0.786610186100006
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.20243358612061
    Y = -0.873619139194489
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = -0.95105654001236
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.40125858783722
    Y = -1.01807391643524
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = -1.07393705844879
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.53884172439575
    Y = -1.1180340051651
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58267605304718
    Y = -1.14988148212433
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.60917019844055
    Y = -1.16913056373596
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6180340051651
    Y = -1.17557048797607
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.60917019844055
    Y = -1.16913056373596
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58267605304718
    Y = -1.14988148212433
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.53884172439575
    Y = -1.1180340051651
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = -1.07393705844879
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.40125858783722
    Y = -1.01807391643524
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = -0.95105654001236
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.20243358612061
    Y = -0.873619139194489
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.08267605304718
    Y = -0.786610186100006
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = -0.69098299741745
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = -0.587785243988037
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.658113718032837
    Y = -0.478147596120834
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = -0.363271266222
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.336408168077469
    Y = -0.24441485106945
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = -0.122880578041077
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.1553595662117
    Y = 0.139886394143105
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = 0.278240144252777
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.459288746118546
    Y = 0.413545459508896
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = 0.544319868087769
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = 0.669130623340607
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.873619139194489
    Y = 0.786610186100006
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = 0.895471513271332
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.10452842712402
    Y = 0.994521915912628
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.20243358612061
    Y = 1.08267605304718
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.28716456890106
    Y = 1.15896821022034
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.35779321193695
    Y = 1.22256243228912
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.41354548931122
    Y = 1.27276206016541
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.4538106918335
    Y = 1.30901694297791
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = 1.33093011379242
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.48628962039948
    Y = 1.33826124668121
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = 1.33093011379242
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.4538106918335
    Y = 1.30901694297791
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.41354548931122
    Y = 1.27276206016541
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.35779321193695
    Y = 1.22256243228912
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.28716456890106
    Y = 1.15896821022034
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.20243358612061
    Y = 1.08267605304718
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.10452842712402
    Y = 0.994521915912628
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = 0.895471513271332
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.873619139194489
    Y = 0.786610186100006
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = 0.669130623340607
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = 0.544319868087769
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.459288746118546
    Y = 0.413545459508896
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = 0.278240144252777
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.1553595662117
    Y = 0.139886394143105
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.82017983593191E-16
    Y = 1.63889731516495E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.1553595662117
    Y = -0.139886394143105
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = -0.278240144252777
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.459288746118546
    Y = -0.413545459508896
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = -0.544319868087769
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = -0.669130623340607
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.873619139194489
    Y = -0.786610186100006
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = -0.895471513271332
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.10452842712402
    Y = -0.994521915912628
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.20243358612061
    Y = -1.08267605304718
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.28716456890106
    Y = -1.15896821022034
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.35779321193695
    Y = -1.22256243228912
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.41354548931122
    Y = -1.27276206016541
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.4538106918335
    Y = -1.30901694297791
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = -1.33093011379242
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.48628962039948
    Y = -1.33826124668121
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = -1.33093011379242
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.4538106918335
    Y = -1.30901694297791
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.41354548931122
    Y = -1.27276206016541
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.35779321193695
    Y = -1.22256243228912
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.28716456890106
    Y = -1.15896821022034
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.20243358612061
    Y = -1.08267605304718
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.10452842712402
    Y = -0.994521915912628
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = -0.895471513271332
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.873619139194489
    Y = -0.786610186100006
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = -0.669130623340607
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = -0.544319868087769
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.459288746118546
    Y = -0.413545459508896
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = -0.278240144252777
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.1553595662117
    Y = -0.139886394143105
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.139886394143105
    Y = 0.1553595662117
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.278240144252777
    Y = 0.30901700258255
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = 0.459288746118546
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.544319868087769
    Y = 0.604528486728668
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.669130623340607
    Y = 0.743144810199738
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.786610186100006
    Y = 0.873619139194489
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.895471513271332
    Y = 0.994521915912628
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = 1.10452842712402
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.08267605304718
    Y = 1.20243358612061
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.15896821022034
    Y = 1.28716456890106
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.22256243228912
    Y = 1.35779321193695
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.27276206016541
    Y = 1.41354548931122
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = 1.4538106918335
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33093011379242
    Y = 1.47814762592316
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33826124668121
    Y = 1.48628962039948
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33093011379242
    Y = 1.47814762592316
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = 1.4538106918335
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.27276206016541
    Y = 1.41354548931122
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.22256243228912
    Y = 1.35779321193695
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.15896821022034
    Y = 1.28716456890106
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.08267605304718
    Y = 1.20243358612061
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = 1.10452842712402
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.895471513271332
    Y = 0.994521915912628
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.786610186100006
    Y = 0.873619139194489
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.669130623340607
    Y = 0.743144810199738
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.544319868087769
    Y = 0.604528486728668
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = 0.459288746118546
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.278240144252777
    Y = 0.30901700258255
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.139886394143105
    Y = 0.1553595662117
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.63889731516495E-16
    Y = 1.82017983593191E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.139886394143105
    Y = -0.1553595662117
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.278240144252777
    Y = -0.30901700258255
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = -0.459288746118546
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.544319868087769
    Y = -0.604528486728668
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.669130623340607
    Y = -0.743144810199738
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.786610186100006
    Y = -0.873619139194489
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.895471513271332
    Y = -0.994521915912628
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = -1.10452842712402
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.08267605304718
    Y = -1.20243358612061
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.15896821022034
    Y = -1.28716456890106
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.22256243228912
    Y = -1.35779321193695
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.27276206016541
    Y = -1.41354548931122
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = -1.4538106918335
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33093011379242
    Y = -1.47814762592316
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33826124668121
    Y = -1.48628962039948
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33093011379242
    Y = -1.47814762592316
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = -1.4538106918335
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.27276206016541
    Y = -1.41354548931122
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.22256243228912
    Y = -1.35779321193695
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.15896821022034
    Y = -1.28716456890106
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.08267605304718
    Y = -1.20243358612061
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = -1.10452842712402
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.895471513271332
    Y = -0.994521915912628
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.786610186100006
    Y = -0.873619139194489
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.669130623340607
    Y = -0.743144810199738
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.544319868087769
    Y = -0.604528486728668
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = -0.459288746118546
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.278240144252777
    Y = -0.30901700258255
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.139886394143105
    Y = -0.1553595662117
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.122880578041077
    Y = 0.169130608439445
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.24441485106945
    Y = 0.336408168077469
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.363271266222
    Y = 0.5
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.478147596120834
    Y = 0.658113718032837
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = 0.80901700258255
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.69098299741745
    Y = 0.95105654001236
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.786610186100006
    Y = 1.08267605304718
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.873619139194489
    Y = 1.20243358612061
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = 1.30901694297791
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.01807391643524
    Y = 1.40125858783722
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.07393705844879
    Y = 1.47814762592316
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.1180340051651
    Y = 1.53884172439575
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.14988148212433
    Y = 1.58267605304718
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.16913056373596
    Y = 1.60917019844055
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.17557048797607
    Y = 1.6180340051651
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.16913056373596
    Y = 1.60917019844055
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.14988148212433
    Y = 1.58267605304718
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.1180340051651
    Y = 1.53884172439575
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.07393705844879
    Y = 1.47814762592316
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.01807391643524
    Y = 1.40125858783722
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = 1.30901694297791
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.873619139194489
    Y = 1.20243358612061
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.786610186100006
    Y = 1.08267605304718
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.69098299741745
    Y = 0.95105654001236
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = 0.80901700258255
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.478147596120834
    Y = 0.658113718032837
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.363271266222
    Y = 0.5
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.24441485106945
    Y = 0.336408168077469
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.122880578041077
    Y = 0.169130608439445
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.43965862235849E-16
    Y = 1.98152015709584E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.122880578041077
    Y = -0.169130608439445
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.24441485106945
    Y = -0.336408168077469
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.363271266222
    Y = -0.5
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.478147596120834
    Y = -0.658113718032837
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = -0.80901700258255
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.69098299741745
    Y = -0.95105654001236
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.786610186100006
    Y = -1.08267605304718
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.873619139194489
    Y = -1.20243358612061
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = -1.30901694297791
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.01807391643524
    Y = -1.40125858783722
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.07393705844879
    Y = -1.47814762592316
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.1180340051651
    Y = -1.53884172439575
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.14988148212433
    Y = -1.58267605304718
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.16913056373596
    Y = -1.60917019844055
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.17557048797607
    Y = -1.6180340051651
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.16913056373596
    Y = -1.60917019844055
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.14988148212433
    Y = -1.58267605304718
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.1180340051651
    Y = -1.53884172439575
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.07393705844879
    Y = -1.47814762592316
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.01807391643524
    Y = -1.40125858783722
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = -1.30901694297791
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.873619139194489
    Y = -1.20243358612061
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.786610186100006
    Y = -1.08267605304718
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.69098299741745
    Y = -0.95105654001236
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = -0.80901700258255
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.478147596120834
    Y = -0.658113718032837
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.363271266222
    Y = -0.5
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.24441485106945
    Y = -0.336408168077469
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.122880578041077
    Y = -0.169130608439445
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.104528464376926
    Y = 0.181048601865768
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = 0.360113620758057
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = 0.53523313999176
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = 0.7044885158538819
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = 0.866025388240814
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = 1.01807391643524
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.669130623340607
    Y = 1.15896821022034
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = 1.28716456890106
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = 1.40125858783722
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.866025388240814
    Y = 1.5
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.913545429706573
    Y = 1.58230710029602
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = 1.64727818965912
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.978147625923157
    Y = 1.6942013502121
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = 1.72256243228912
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1
    Y = 1.73205077648163
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = 1.72256243228912
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.978147625923157
    Y = 1.6942013502121
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = 1.64727818965912
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.913545429706573
    Y = 1.58230710029602
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.866025388240814
    Y = 1.5
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = 1.40125858783722
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = 1.28716456890106
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.669130623340607
    Y = 1.15896821022034
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = 1.01807391643524
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = 0.866025388240814
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = 0.7044885158538819
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = 0.53523313999176
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = 0.360113620758057
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.104528464376926
    Y = 0.181048601865768
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.22464685258517E-16
    Y = 2.12115049442606E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.104528464376926
    Y = -0.181048601865768
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = -0.360113620758057
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = -0.53523313999176
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = -0.7044885158538819
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = -0.866025388240814
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = -1.01807391643524
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.669130623340607
    Y = -1.15896821022034
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = -1.28716456890106
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = -1.40125858783722
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.866025388240814
    Y = -1.5
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.913545429706573
    Y = -1.58230710029602
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = -1.64727818965912
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.978147625923157
    Y = -1.6942013502121
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = -1.72256243228912
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1
    Y = -1.73205077648163
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = -1.72256243228912
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.978147625923157
    Y = -1.6942013502121
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = -1.64727818965912
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.913545429706573
    Y = -1.58230710029602
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.866025388240814
    Y = -1.5
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = -1.40125858783722
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = -1.28716456890106
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.669130623340607
    Y = -1.15896821022034
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = -1.01807391643524
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = -0.866025388240814
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = -0.7044885158538819
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = -0.53523313999176
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = -0.360113620758057
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.104528464376926
    Y = -0.181048601865768
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0850311145186424
    Y = 0.190983012318611
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = 0.37987357378006
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.251377075910568
    Y = 0.564602136611938
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.330869406461716
    Y = 0.743144810199738
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = 0.913545429706573
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.478147596120834
    Y = 1.07393705844879
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.544319868087769
    Y = 1.22256243228912
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = 1.35779321193695
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.658113718032837
    Y = 1.47814762592316
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.7044885158538819
    Y = 1.58230710029602
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = 1.66913056373596
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.77365905046463
    Y = 1.73766672611237
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.795696914196014
    Y = 1.78716456890106
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = 1.81708192825317
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.813473284244537
    Y = 1.82709085941315
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = 1.81708192825317
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.795696914196014
    Y = 1.78716456890106
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.77365905046463
    Y = 1.73766672611237
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = 1.66913056373596
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.7044885158538819
    Y = 1.58230710029602
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.658113718032837
    Y = 1.47814762592316
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = 1.35779321193695
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.544319868087769
    Y = 1.22256243228912
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.478147596120834
    Y = 1.07393705844879
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = 0.913545429706573
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.330869406461716
    Y = 0.743144810199738
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.251377075910568
    Y = 0.564602136611938
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = 0.37987357378006
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0850311145186424
    Y = 0.190983012318611
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 9.96217485357301E-17
    Y = 2.23754102701048E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0850311145186424
    Y = -0.190983012318611
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = -0.37987357378006
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.251377075910568
    Y = -0.564602136611938
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.330869406461716
    Y = -0.743144810199738
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = -0.913545429706573
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.478147596120834
    Y = -1.07393705844879
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.544319868087769
    Y = -1.22256243228912
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = -1.35779321193695
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.658113718032837
    Y = -1.47814762592316
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.7044885158538819
    Y = -1.58230710029602
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = -1.66913056373596
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.77365905046463
    Y = -1.73766672611237
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.795696914196014
    Y = -1.78716456890106
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = -1.81708192825317
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.813473284244537
    Y = -1.82709085941315
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = -1.81708192825317
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.795696914196014
    Y = -1.78716456890106
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.77365905046463
    Y = -1.73766672611237
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = -1.66913056373596
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.7044885158538819
    Y = -1.58230710029602
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.658113718032837
    Y = -1.47814762592316
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = -1.35779321193695
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.544319868087769
    Y = -1.22256243228912
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.478147596120834
    Y = -1.07393705844879
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = -0.913545429706573
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.330869406461716
    Y = -0.743144810199738
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.251377075910568
    Y = -0.564602136611938
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = -0.37987357378006
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0850311145186424
    Y = -0.190983012318611
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0646021440625191
    Y = 0.19882495701313
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.128496497869492
    Y = 0.395471543073654
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = 0.587785243988037
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.251377075910568
    Y = 0.77365905046463
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = 0.95105654001236
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.363271266222
    Y = 1.1180340051651
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = 1.27276206016541
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.459288746118546
    Y = 1.41354548931122
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = 1.53884172439575
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.53523313999176
    Y = 1.64727818965912
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.564602136611938
    Y = 1.73766672611237
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = 1.80901694297791
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = 1.86054730415344
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.614648342132568
    Y = 1.89169311523438
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.6180340051651
    Y = 1.90211308002472
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.614648342132568
    Y = 1.89169311523438
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = 1.86054730415344
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = 1.80901694297791
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.564602136611938
    Y = 1.73766672611237
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.53523313999176
    Y = 1.64727818965912
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = 1.53884172439575
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.459288746118546
    Y = 1.41354548931122
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = 1.27276206016541
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.363271266222
    Y = 1.1180340051651
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = 0.95105654001236
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.251377075910568
    Y = 0.77365905046463
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = 0.587785243988037
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.128496497869492
    Y = 0.395471543073654
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0646021440625191
    Y = 0.19882495701313
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 7.56873370685124E-17
    Y = 2.32941657321678E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0646021440625191
    Y = -0.19882495701313
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.128496497869492
    Y = -0.395471543073654
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = -0.587785243988037
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.251377075910568
    Y = -0.77365905046463
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = -0.95105654001236
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.363271266222
    Y = -1.1180340051651
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = -1.27276206016541
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.459288746118546
    Y = -1.41354548931122
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = -1.53884172439575
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.53523313999176
    Y = -1.64727818965912
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.564602136611938
    Y = -1.73766672611237
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = -1.80901694297791
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = -1.86054730415344
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.614648342132568
    Y = -1.89169311523438
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.6180340051651
    Y = -1.90211308002472
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.614648342132568
    Y = -1.89169311523438
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = -1.86054730415344
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = -1.80901694297791
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.564602136611938
    Y = -1.73766672611237
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.53523313999176
    Y = -1.64727818965912
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = -1.53884172439575
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.459288746118546
    Y = -1.41354548931122
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = -1.27276206016541
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.363271266222
    Y = -1.1180340051651
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = -0.95105654001236
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.251377075910568
    Y = -0.77365905046463
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = -0.587785243988037
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.128496497869492
    Y = -0.395471543073654
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0646021440625191
    Y = -0.19882495701313
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0434653796255589
    Y = 0.204488530755043
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0864545404911041
    Y = 0.406736642122269
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.128496497869492
    Y = 0.604528486728668
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = 0.795696914196014
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = 0.978147625923157
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.24441485106945
    Y = 1.14988148212433
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.278240144252777
    Y = 1.30901694297791
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = 1.4538106918335
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.336408168077469
    Y = 1.58267605304718
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.360113620758057
    Y = 1.6942013502121
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.37987357378006
    Y = 1.78716456890106
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.395471543073654
    Y = 1.86054730415344
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = 1.91354548931122
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = 1.94557845592499
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.415823370218277
    Y = 1.95629525184631
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = 1.94557845592499
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = 1.91354548931122
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.395471543073654
    Y = 1.86054730415344
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.37987357378006
    Y = 1.78716456890106
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.360113620758057
    Y = 1.6942013502121
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.336408168077469
    Y = 1.58267605304718
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = 1.4538106918335
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.278240144252777
    Y = 1.30901694297791
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.24441485106945
    Y = 1.14988148212433
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = 0.978147625923157
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = 0.795696914196014
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.128496497869492
    Y = 0.604528486728668
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0864545404911041
    Y = 0.406736642122269
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0434653796255589
    Y = 0.204488530755043
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 5.09236770372102E-17
    Y = 2.39577075202452E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0434653796255589
    Y = -0.204488530755043
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0864545404911041
    Y = -0.406736642122269
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.128496497869492
    Y = -0.604528486728668
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = -0.795696914196014
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = -0.978147625923157
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.24441485106945
    Y = -1.14988148212433
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.278240144252777
    Y = -1.30901694297791
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = -1.4538106918335
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.336408168077469
    Y = -1.58267605304718
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.360113620758057
    Y = -1.6942013502121
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.37987357378006
    Y = -1.78716456890106
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.395471543073654
    Y = -1.86054730415344
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = -1.91354548931122
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = -1.94557845592499
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.415823370218277
    Y = -1.95629525184631
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = -1.94557845592499
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = -1.91354548931122
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.395471543073654
    Y = -1.86054730415344
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.37987357378006
    Y = -1.78716456890106
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.360113620758057
    Y = -1.6942013502121
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.336408168077469
    Y = -1.58267605304718
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = -1.4538106918335
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.278240144252777
    Y = -1.30901694297791
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.24441485106945
    Y = -1.14988148212433
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = -0.978147625923157
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = -0.795696914196014
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.128496497869492
    Y = -0.604528486728668
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0864545404911041
    Y = -0.406736642122269
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0434653796255589
    Y = -0.204488530755043
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0218524001538754
    Y = 0.207911685109138
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0434653796255589
    Y = 0.413545459508896
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0646021440625191
    Y = 0.614648342132568
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0850311145186424
    Y = 0.80901700258255
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.104528464376926
    Y = 0.994521915912628
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.122880578041077
    Y = 1.16913056373596
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.139886394143105
    Y = 1.33093011379242
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.1553595662117
    Y = 1.47814762592316
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = 1.60917019844055
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.181048601865768
    Y = 1.72256243228912
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = 1.81708192825317
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.19882495701313
    Y = 1.89169311523438
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.204488530755043
    Y = 1.94557845592499
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = 1.97814762592316
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.209056928753853
    Y = 1.98904383182526
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = 1.97814762592316
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.204488530755043
    Y = 1.94557845592499
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.19882495701313
    Y = 1.89169311523438
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = 1.81708192825317
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.181048601865768
    Y = 1.72256243228912
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = 1.60917019844055
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.1553595662117
    Y = 1.47814762592316
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.139886394143105
    Y = 1.33093011379242
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.122880578041077
    Y = 1.16913056373596
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.104528464376926
    Y = 0.994521915912628
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0850311145186424
    Y = 0.80901700258255
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0646021440625191
    Y = 0.614648342132568
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0434653796255589
    Y = 0.413545459508896
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0218524001538754
    Y = 0.207911685109138
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 2.56020903001859E-17
    Y = 2.43587617389024E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0218524001538754
    Y = -0.207911685109138
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0434653796255589
    Y = -0.413545459508896
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0646021440625191
    Y = -0.614648342132568
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0850311145186424
    Y = -0.80901700258255
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.104528464376926
    Y = -0.994521915912628
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.122880578041077
    Y = -1.16913056373596
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.139886394143105
    Y = -1.33093011379242
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.1553595662117
    Y = -1.47814762592316
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = -1.60917019844055
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.181048601865768
    Y = -1.72256243228912
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = -1.81708192825317
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.19882495701313
    Y = -1.89169311523438
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.204488530755043
    Y = -1.94557845592499
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = -1.97814762592316
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.209056928753853
    Y = -1.98904383182526
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = -1.97814762592316
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.204488530755043
    Y = -1.94557845592499
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.19882495701313
    Y = -1.89169311523438
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = -1.81708192825317
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.181048601865768
    Y = -1.72256243228912
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = -1.60917019844055
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.1553595662117
    Y = -1.47814762592316
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.139886394143105
    Y = -1.33093011379242
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.122880578041077
    Y = -1.16913056373596
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.104528464376926
    Y = -0.994521915912628
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0850311145186424
    Y = -0.80901700258255
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0646021440625191
    Y = -0.614648342132568
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0434653796255589
    Y = -0.413545459508896
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0218524001538754
    Y = -0.207911685109138
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.2801045150093E-17
    Y = 0.209056928753853
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 2.54618385186051E-17
    Y = 0.415823370218277
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 3.78436685342562E-17
    Y = 0.6180340051651
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 4.98108742678651E-17
    Y = 0.813473284244537
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 6.123234262925841E-17
    Y = 1
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 7.19829311179243E-17
    Y = 1.17557048797607
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 8.19448657582474E-17
    Y = 1.33826124668121
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 9.100899179659529E-17
    Y = 1.48628962039948
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 9.907600785479209E-17
    Y = 1.6180340051651
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.06057524721303E-16
    Y = 1.73205077648163
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.11877051350524E-16
    Y = 1.82709085941315
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.16470828660839E-16
    Y = 1.90211308002472
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.19788537601226E-16
    Y = 1.95629525184631
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.21793808694512E-16
    Y = 1.98904383182526
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.22464685258517E-16
    Y = 2
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.21793808694512E-16
    Y = 1.98904383182526
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.19788537601226E-16
    Y = 1.95629525184631
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.16470828660839E-16
    Y = 1.90211308002472
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.11877051350524E-16
    Y = 1.82709085941315
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.06057524721303E-16
    Y = 1.73205077648163
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 9.907600785479209E-17
    Y = 1.6180340051651
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 9.100899179659529E-17
    Y = 1.48628962039948
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 8.19448657582474E-17
    Y = 1.33826124668121
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 7.19829311179243E-17
    Y = 1.17557048797607
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 6.123234262925841E-17
    Y = 1
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 4.98108742678651E-17
    Y = 0.813473284244537
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 3.78436685342562E-17
    Y = 0.6180340051651
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 2.54618385186051E-17
    Y = 0.415823370218277
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.2801045150093E-17
    Y = 0.209056928753853
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.49975975722119E-32
    Y = 2.44929370517034E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.2801045150093E-17
    Y = -0.209056928753853
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -2.54618385186051E-17
    Y = -0.415823370218277
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -3.78436685342562E-17
    Y = -0.6180340051651
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -4.98108742678651E-17
    Y = -0.813473284244537
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -6.123234262925841E-17
    Y = -1
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -7.19829311179243E-17
    Y = -1.17557048797607
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -8.19448657582474E-17
    Y = -1.33826124668121
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -9.100899179659529E-17
    Y = -1.48628962039948
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -9.907600785479209E-17
    Y = -1.6180340051651
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.06057524721303E-16
    Y = -1.73205077648163
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.11877051350524E-16
    Y = -1.82709085941315
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.16470828660839E-16
    Y = -1.90211308002472
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.19788537601226E-16
    Y = -1.95629525184631
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.21793808694512E-16
    Y = -1.98904383182526
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.22464685258517E-16
    Y = -2
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.21793808694512E-16
    Y = -1.98904383182526
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.19788537601226E-16
    Y = -1.95629525184631
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.16470828660839E-16
    Y = -1.90211308002472
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.11877051350524E-16
    Y = -1.82709085941315
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.06057524721303E-16
    Y = -1.73205077648163
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -9.907600785479209E-17
    Y = -1.6180340051651
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -9.100899179659529E-17
    Y = -1.48628962039948
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -8.19448657582474E-17
    Y = -1.33826124668121
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -7.19829311179243E-17
    Y = -1.17557048797607
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -6.123234262925841E-17
    Y = -1
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -4.98108742678651E-17
    Y = -0.813473284244537
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -3.78436685342562E-17
    Y = -0.6180340051651
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -2.54618385186051E-17
    Y = -0.415823370218277
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.2801045150093E-17
    Y = -0.209056928753853
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0218524001538754
    Y = 0.207911685109138
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0434653796255589
    Y = 0.413545459508896
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0646021440625191
    Y = 0.614648342132568
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0850311145186424
    Y = 0.80901700258255
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.104528464376926
    Y = 0.994521915912628
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.122880578041077
    Y = 1.16913056373596
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.139886394143105
    Y = 1.33093011379242
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.1553595662117
    Y = 1.47814762592316
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = 1.60917019844055
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.181048601865768
    Y = 1.72256243228912
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = 1.81708192825317
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.19882495701313
    Y = 1.89169311523438
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.204488530755043
    Y = 1.94557845592499
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = 1.97814762592316
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.209056928753853
    Y = 1.98904383182526
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = 1.97814762592316
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.204488530755043
    Y = 1.94557845592499
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.19882495701313
    Y = 1.89169311523438
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = 1.81708192825317
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.181048601865768
    Y = 1.72256243228912
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = 1.60917019844055
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.1553595662117
    Y = 1.47814762592316
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.139886394143105
    Y = 1.33093011379242
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.122880578041077
    Y = 1.16913056373596
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.104528464376926
    Y = 0.994521915912628
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0850311145186424
    Y = 0.80901700258255
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0646021440625191
    Y = 0.614648342132568
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0434653796255589
    Y = 0.413545459508896
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0218524001538754
    Y = 0.207911685109138
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -2.56020903001859E-17
    Y = 2.43587617389024E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0218524001538754
    Y = -0.207911685109138
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0434653796255589
    Y = -0.413545459508896
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0646021440625191
    Y = -0.614648342132568
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0850311145186424
    Y = -0.80901700258255
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.104528464376926
    Y = -0.994521915912628
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.122880578041077
    Y = -1.16913056373596
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.139886394143105
    Y = -1.33093011379242
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.1553595662117
    Y = -1.47814762592316
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = -1.60917019844055
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.181048601865768
    Y = -1.72256243228912
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = -1.81708192825317
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.19882495701313
    Y = -1.89169311523438
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.204488530755043
    Y = -1.94557845592499
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = -1.97814762592316
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.209056928753853
    Y = -1.98904383182526
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = -1.97814762592316
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.204488530755043
    Y = -1.94557845592499
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.19882495701313
    Y = -1.89169311523438
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = -1.81708192825317
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.181048601865768
    Y = -1.72256243228912
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = -1.60917019844055
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.1553595662117
    Y = -1.47814762592316
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.139886394143105
    Y = -1.33093011379242
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.122880578041077
    Y = -1.16913056373596
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.104528464376926
    Y = -0.994521915912628
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0850311145186424
    Y = -0.80901700258255
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0646021440625191
    Y = -0.614648342132568
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0434653796255589
    Y = -0.413545459508896
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0218524001538754
    Y = -0.207911685109138
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0434653796255589
    Y = 0.204488530755043
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0864545404911041
    Y = 0.406736642122269
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.128496497869492
    Y = 0.604528486728668
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = 0.795696914196014
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = 0.978147625923157
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.24441485106945
    Y = 1.14988148212433
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.278240144252777
    Y = 1.30901694297791
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = 1.4538106918335
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.336408168077469
    Y = 1.58267605304718
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.360113620758057
    Y = 1.6942013502121
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.37987357378006
    Y = 1.78716456890106
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.395471543073654
    Y = 1.86054730415344
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = 1.91354548931122
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = 1.94557845592499
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.415823370218277
    Y = 1.95629525184631
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = 1.94557845592499
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = 1.91354548931122
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.395471543073654
    Y = 1.86054730415344
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.37987357378006
    Y = 1.78716456890106
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.360113620758057
    Y = 1.6942013502121
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.336408168077469
    Y = 1.58267605304718
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = 1.4538106918335
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.278240144252777
    Y = 1.30901694297791
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.24441485106945
    Y = 1.14988148212433
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = 0.978147625923157
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = 0.795696914196014
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.128496497869492
    Y = 0.604528486728668
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0864545404911041
    Y = 0.406736642122269
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0434653796255589
    Y = 0.204488530755043
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -5.09236770372102E-17
    Y = 2.39577075202452E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0434653796255589
    Y = -0.204488530755043
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0864545404911041
    Y = -0.406736642122269
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.128496497869492
    Y = -0.604528486728668
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = -0.795696914196014
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = -0.978147625923157
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.24441485106945
    Y = -1.14988148212433
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.278240144252777
    Y = -1.30901694297791
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = -1.4538106918335
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.336408168077469
    Y = -1.58267605304718
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.360113620758057
    Y = -1.6942013502121
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.37987357378006
    Y = -1.78716456890106
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.395471543073654
    Y = -1.86054730415344
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = -1.91354548931122
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = -1.94557845592499
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.415823370218277
    Y = -1.95629525184631
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = -1.94557845592499
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = -1.91354548931122
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.395471543073654
    Y = -1.86054730415344
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.37987357378006
    Y = -1.78716456890106
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.360113620758057
    Y = -1.6942013502121
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.336408168077469
    Y = -1.58267605304718
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = -1.4538106918335
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.278240144252777
    Y = -1.30901694297791
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.24441485106945
    Y = -1.14988148212433
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = -0.978147625923157
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = -0.795696914196014
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.128496497869492
    Y = -0.604528486728668
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0864545404911041
    Y = -0.406736642122269
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0434653796255589
    Y = -0.204488530755043
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0646021440625191
    Y = 0.19882495701313
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.128496497869492
    Y = 0.395471543073654
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = 0.587785243988037
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.251377075910568
    Y = 0.77365905046463
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = 0.95105654001236
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.363271266222
    Y = 1.1180340051651
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = 1.27276206016541
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.459288746118546
    Y = 1.41354548931122
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = 1.53884172439575
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.53523313999176
    Y = 1.64727818965912
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.564602136611938
    Y = 1.73766672611237
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = 1.80901694297791
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = 1.86054730415344
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.614648342132568
    Y = 1.89169311523438
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.6180340051651
    Y = 1.90211308002472
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.614648342132568
    Y = 1.89169311523438
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = 1.86054730415344
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = 1.80901694297791
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.564602136611938
    Y = 1.73766672611237
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.53523313999176
    Y = 1.64727818965912
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = 1.53884172439575
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.459288746118546
    Y = 1.41354548931122
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = 1.27276206016541
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.363271266222
    Y = 1.1180340051651
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = 0.95105654001236
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.251377075910568
    Y = 0.77365905046463
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = 0.587785243988037
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.128496497869492
    Y = 0.395471543073654
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0646021440625191
    Y = 0.19882495701313
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -7.56873370685124E-17
    Y = 2.32941657321678E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0646021440625191
    Y = -0.19882495701313
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.128496497869492
    Y = -0.395471543073654
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = -0.587785243988037
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.251377075910568
    Y = -0.77365905046463
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = -0.95105654001236
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.363271266222
    Y = -1.1180340051651
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = -1.27276206016541
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.459288746118546
    Y = -1.41354548931122
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = -1.53884172439575
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.53523313999176
    Y = -1.64727818965912
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.564602136611938
    Y = -1.73766672611237
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = -1.80901694297791
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = -1.86054730415344
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.614648342132568
    Y = -1.89169311523438
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.6180340051651
    Y = -1.90211308002472
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.614648342132568
    Y = -1.89169311523438
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = -1.86054730415344
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = -1.80901694297791
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.564602136611938
    Y = -1.73766672611237
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.53523313999176
    Y = -1.64727818965912
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = -1.53884172439575
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.459288746118546
    Y = -1.41354548931122
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = -1.27276206016541
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.363271266222
    Y = -1.1180340051651
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = -0.95105654001236
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.251377075910568
    Y = -0.77365905046463
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = -0.587785243988037
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.128496497869492
    Y = -0.395471543073654
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0646021440625191
    Y = -0.19882495701313
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0850311145186424
    Y = 0.190983012318611
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = 0.37987357378006
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.251377075910568
    Y = 0.564602136611938
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.330869406461716
    Y = 0.743144810199738
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = 0.913545429706573
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.478147596120834
    Y = 1.07393705844879
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.544319868087769
    Y = 1.22256243228912
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = 1.35779321193695
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.658113718032837
    Y = 1.47814762592316
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.7044885158538819
    Y = 1.58230710029602
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = 1.66913056373596
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.77365905046463
    Y = 1.73766672611237
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.795696914196014
    Y = 1.78716456890106
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = 1.81708192825317
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.813473284244537
    Y = 1.82709085941315
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = 1.81708192825317
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.795696914196014
    Y = 1.78716456890106
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.77365905046463
    Y = 1.73766672611237
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = 1.66913056373596
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.7044885158538819
    Y = 1.58230710029602
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.658113718032837
    Y = 1.47814762592316
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = 1.35779321193695
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.544319868087769
    Y = 1.22256243228912
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.478147596120834
    Y = 1.07393705844879
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = 0.913545429706573
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.330869406461716
    Y = 0.743144810199738
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.251377075910568
    Y = 0.564602136611938
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = 0.37987357378006
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.0850311145186424
    Y = 0.190983012318611
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -9.96217485357301E-17
    Y = 2.23754102701048E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0850311145186424
    Y = -0.190983012318611
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = -0.37987357378006
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.251377075910568
    Y = -0.564602136611938
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.330869406461716
    Y = -0.743144810199738
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = -0.913545429706573
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.478147596120834
    Y = -1.07393705844879
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.544319868087769
    Y = -1.22256243228912
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = -1.35779321193695
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.658113718032837
    Y = -1.47814762592316
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.7044885158538819
    Y = -1.58230710029602
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = -1.66913056373596
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.77365905046463
    Y = -1.73766672611237
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.795696914196014
    Y = -1.78716456890106
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = -1.81708192825317
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.813473284244537
    Y = -1.82709085941315
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = -1.81708192825317
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.795696914196014
    Y = -1.78716456890106
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.77365905046463
    Y = -1.73766672611237
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = -1.66913056373596
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.7044885158538819
    Y = -1.58230710029602
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.658113718032837
    Y = -1.47814762592316
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = -1.35779321193695
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.544319868087769
    Y = -1.22256243228912
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.478147596120834
    Y = -1.07393705844879
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = -0.913545429706573
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.330869406461716
    Y = -0.743144810199738
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.251377075910568
    Y = -0.564602136611938
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = -0.37987357378006
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.0850311145186424
    Y = -0.190983012318611
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.104528464376926
    Y = 0.181048601865768
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = 0.360113620758057
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = 0.53523313999176
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = 0.7044885158538819
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = 0.866025388240814
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = 1.01807391643524
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.669130623340607
    Y = 1.15896821022034
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = 1.28716456890106
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = 1.40125858783722
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.866025388240814
    Y = 1.5
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.913545429706573
    Y = 1.58230710029602
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = 1.64727818965912
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.978147625923157
    Y = 1.6942013502121
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = 1.72256243228912
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1
    Y = 1.73205077648163
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = 1.72256243228912
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.978147625923157
    Y = 1.6942013502121
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = 1.64727818965912
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.913545429706573
    Y = 1.58230710029602
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.866025388240814
    Y = 1.5
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = 1.40125858783722
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = 1.28716456890106
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.669130623340607
    Y = 1.15896821022034
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = 1.01807391643524
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = 0.866025388240814
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = 0.7044885158538819
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = 0.53523313999176
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = 0.360113620758057
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.104528464376926
    Y = 0.181048601865768
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.22464685258517E-16
    Y = 2.12115049442606E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.104528464376926
    Y = -0.181048601865768
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = -0.360113620758057
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = -0.53523313999176
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = -0.7044885158538819
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = -0.866025388240814
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = -1.01807391643524
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.669130623340607
    Y = -1.15896821022034
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = -1.28716456890106
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = -1.40125858783722
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.866025388240814
    Y = -1.5
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.913545429706573
    Y = -1.58230710029602
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = -1.64727818965912
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.978147625923157
    Y = -1.6942013502121
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = -1.72256243228912
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1
    Y = -1.73205077648163
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = -1.72256243228912
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.978147625923157
    Y = -1.6942013502121
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = -1.64727818965912
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.913545429706573
    Y = -1.58230710029602
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.866025388240814
    Y = -1.5
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = -1.40125858783722
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = -1.28716456890106
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.669130623340607
    Y = -1.15896821022034
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = -1.01807391643524
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = -0.866025388240814
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = -0.7044885158538819
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = -0.53523313999176
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = -0.360113620758057
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.104528464376926
    Y = -0.181048601865768
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.122880578041077
    Y = 0.169130608439445
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.24441485106945
    Y = 0.336408168077469
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.363271266222
    Y = 0.5
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.478147596120834
    Y = 0.658113718032837
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = 0.80901700258255
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.69098299741745
    Y = 0.95105654001236
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.786610186100006
    Y = 1.08267605304718
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.873619139194489
    Y = 1.20243358612061
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = 1.30901694297791
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.01807391643524
    Y = 1.40125858783722
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.07393705844879
    Y = 1.47814762592316
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.1180340051651
    Y = 1.53884172439575
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.14988148212433
    Y = 1.58267605304718
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.16913056373596
    Y = 1.60917019844055
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.17557048797607
    Y = 1.6180340051651
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.16913056373596
    Y = 1.60917019844055
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.14988148212433
    Y = 1.58267605304718
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.1180340051651
    Y = 1.53884172439575
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.07393705844879
    Y = 1.47814762592316
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.01807391643524
    Y = 1.40125858783722
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = 1.30901694297791
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.873619139194489
    Y = 1.20243358612061
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.786610186100006
    Y = 1.08267605304718
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.69098299741745
    Y = 0.95105654001236
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = 0.80901700258255
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.478147596120834
    Y = 0.658113718032837
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.363271266222
    Y = 0.5
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.24441485106945
    Y = 0.336408168077469
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.122880578041077
    Y = 0.169130608439445
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.43965862235849E-16
    Y = 1.98152015709584E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.122880578041077
    Y = -0.169130608439445
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.24441485106945
    Y = -0.336408168077469
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.363271266222
    Y = -0.5
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.478147596120834
    Y = -0.658113718032837
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = -0.80901700258255
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.69098299741745
    Y = -0.95105654001236
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.786610186100006
    Y = -1.08267605304718
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.873619139194489
    Y = -1.20243358612061
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = -1.30901694297791
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.01807391643524
    Y = -1.40125858783722
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.07393705844879
    Y = -1.47814762592316
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.1180340051651
    Y = -1.53884172439575
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.14988148212433
    Y = -1.58267605304718
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.16913056373596
    Y = -1.60917019844055
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.17557048797607
    Y = -1.6180340051651
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.16913056373596
    Y = -1.60917019844055
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.14988148212433
    Y = -1.58267605304718
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.1180340051651
    Y = -1.53884172439575
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.07393705844879
    Y = -1.47814762592316
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.01807391643524
    Y = -1.40125858783722
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = -1.30901694297791
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.873619139194489
    Y = -1.20243358612061
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.786610186100006
    Y = -1.08267605304718
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.69098299741745
    Y = -0.95105654001236
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = -0.80901700258255
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.478147596120834
    Y = -0.658113718032837
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.363271266222
    Y = -0.5
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.24441485106945
    Y = -0.336408168077469
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.122880578041077
    Y = -0.169130608439445
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.139886394143105
    Y = 0.1553595662117
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.278240144252777
    Y = 0.30901700258255
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = 0.459288746118546
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.544319868087769
    Y = 0.604528486728668
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.669130623340607
    Y = 0.743144810199738
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.786610186100006
    Y = 0.873619139194489
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.895471513271332
    Y = 0.994521915912628
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = 1.10452842712402
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.08267605304718
    Y = 1.20243358612061
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.15896821022034
    Y = 1.28716456890106
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.22256243228912
    Y = 1.35779321193695
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.27276206016541
    Y = 1.41354548931122
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = 1.4538106918335
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33093011379242
    Y = 1.47814762592316
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33826124668121
    Y = 1.48628962039948
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33093011379242
    Y = 1.47814762592316
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = 1.4538106918335
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.27276206016541
    Y = 1.41354548931122
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.22256243228912
    Y = 1.35779321193695
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.15896821022034
    Y = 1.28716456890106
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.08267605304718
    Y = 1.20243358612061
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = 1.10452842712402
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.895471513271332
    Y = 0.994521915912628
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.786610186100006
    Y = 0.873619139194489
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.669130623340607
    Y = 0.743144810199738
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.544319868087769
    Y = 0.604528486728668
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = 0.459288746118546
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.278240144252777
    Y = 0.30901700258255
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.139886394143105
    Y = 0.1553595662117
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.63889731516495E-16
    Y = 1.82017983593191E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.139886394143105
    Y = -0.1553595662117
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.278240144252777
    Y = -0.30901700258255
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = -0.459288746118546
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.544319868087769
    Y = -0.604528486728668
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.669130623340607
    Y = -0.743144810199738
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.786610186100006
    Y = -0.873619139194489
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.895471513271332
    Y = -0.994521915912628
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = -1.10452842712402
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.08267605304718
    Y = -1.20243358612061
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.15896821022034
    Y = -1.28716456890106
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.22256243228912
    Y = -1.35779321193695
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.27276206016541
    Y = -1.41354548931122
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = -1.4538106918335
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33093011379242
    Y = -1.47814762592316
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33826124668121
    Y = -1.48628962039948
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33093011379242
    Y = -1.47814762592316
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = -1.4538106918335
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.27276206016541
    Y = -1.41354548931122
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.22256243228912
    Y = -1.35779321193695
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.15896821022034
    Y = -1.28716456890106
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.08267605304718
    Y = -1.20243358612061
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = -1.10452842712402
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.895471513271332
    Y = -0.994521915912628
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.786610186100006
    Y = -0.873619139194489
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.669130623340607
    Y = -0.743144810199738
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.544319868087769
    Y = -0.604528486728668
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = -0.459288746118546
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.278240144252777
    Y = -0.30901700258255
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.139886394143105
    Y = -0.1553595662117
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.1553595662117
    Y = 0.139886394143105
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = 0.278240144252777
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.459288746118546
    Y = 0.413545459508896
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = 0.544319868087769
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = 0.669130623340607
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.873619139194489
    Y = 0.786610186100006
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = 0.895471513271332
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.10452842712402
    Y = 0.994521915912628
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.20243358612061
    Y = 1.08267605304718
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.28716456890106
    Y = 1.15896821022034
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.35779321193695
    Y = 1.22256243228912
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.41354548931122
    Y = 1.27276206016541
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.4538106918335
    Y = 1.30901694297791
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = 1.33093011379242
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.48628962039948
    Y = 1.33826124668121
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = 1.33093011379242
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.4538106918335
    Y = 1.30901694297791
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.41354548931122
    Y = 1.27276206016541
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.35779321193695
    Y = 1.22256243228912
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.28716456890106
    Y = 1.15896821022034
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.20243358612061
    Y = 1.08267605304718
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.10452842712402
    Y = 0.994521915912628
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = 0.895471513271332
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.873619139194489
    Y = 0.786610186100006
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = 0.669130623340607
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = 0.544319868087769
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.459288746118546
    Y = 0.413545459508896
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.30901700258255
    Y = 0.278240144252777
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.1553595662117
    Y = 0.139886394143105
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.82017983593191E-16
    Y = 1.63889731516495E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.1553595662117
    Y = -0.139886394143105
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = -0.278240144252777
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.459288746118546
    Y = -0.413545459508896
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = -0.544319868087769
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = -0.669130623340607
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.873619139194489
    Y = -0.786610186100006
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = -0.895471513271332
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.10452842712402
    Y = -0.994521915912628
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.20243358612061
    Y = -1.08267605304718
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.28716456890106
    Y = -1.15896821022034
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.35779321193695
    Y = -1.22256243228912
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.41354548931122
    Y = -1.27276206016541
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.4538106918335
    Y = -1.30901694297791
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = -1.33093011379242
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.48628962039948
    Y = -1.33826124668121
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = -1.33093011379242
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.4538106918335
    Y = -1.30901694297791
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.41354548931122
    Y = -1.27276206016541
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.35779321193695
    Y = -1.22256243228912
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.28716456890106
    Y = -1.15896821022034
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.20243358612061
    Y = -1.08267605304718
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.10452842712402
    Y = -0.994521915912628
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = -0.895471513271332
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.873619139194489
    Y = -0.786610186100006
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = -0.669130623340607
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = -0.544319868087769
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.459288746118546
    Y = -0.413545459508896
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.30901700258255
    Y = -0.278240144252777
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.1553595662117
    Y = -0.139886394143105
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = 0.122880578041077
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.336408168077469
    Y = 0.24441485106945
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = 0.363271266222
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.658113718032837
    Y = 0.478147596120834
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = 0.587785243988037
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = 0.69098299741745
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.08267605304718
    Y = 0.786610186100006
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.20243358612061
    Y = 0.873619139194489
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = 0.95105654001236
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.40125858783722
    Y = 1.01807391643524
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = 1.07393705844879
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.53884172439575
    Y = 1.1180340051651
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58267605304718
    Y = 1.14988148212433
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.60917019844055
    Y = 1.16913056373596
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6180340051651
    Y = 1.17557048797607
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.60917019844055
    Y = 1.16913056373596
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58267605304718
    Y = 1.14988148212433
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.53884172439575
    Y = 1.1180340051651
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = 1.07393705844879
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.40125858783722
    Y = 1.01807391643524
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = 0.95105654001236
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.20243358612061
    Y = 0.873619139194489
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.08267605304718
    Y = 0.786610186100006
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = 0.69098299741745
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = 0.587785243988037
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.658113718032837
    Y = 0.478147596120834
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.5
    Y = 0.363271266222
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.336408168077469
    Y = 0.24441485106945
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.169130608439445
    Y = 0.122880578041077
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.98152015709584E-16
    Y = 1.43965862235849E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = -0.122880578041077
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.336408168077469
    Y = -0.24441485106945
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = -0.363271266222
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.658113718032837
    Y = -0.478147596120834
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = -0.587785243988037
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = -0.69098299741745
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.08267605304718
    Y = -0.786610186100006
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.20243358612061
    Y = -0.873619139194489
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = -0.95105654001236
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.40125858783722
    Y = -1.01807391643524
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = -1.07393705844879
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.53884172439575
    Y = -1.1180340051651
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58267605304718
    Y = -1.14988148212433
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.60917019844055
    Y = -1.16913056373596
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6180340051651
    Y = -1.17557048797607
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.60917019844055
    Y = -1.16913056373596
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58267605304718
    Y = -1.14988148212433
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.53884172439575
    Y = -1.1180340051651
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = -1.07393705844879
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.40125858783722
    Y = -1.01807391643524
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = -0.95105654001236
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.20243358612061
    Y = -0.873619139194489
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.08267605304718
    Y = -0.786610186100006
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = -0.69098299741745
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = -0.587785243988037
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.658113718032837
    Y = -0.478147596120834
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.5
    Y = -0.363271266222
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.336408168077469
    Y = -0.24441485106945
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.169130608439445
    Y = -0.122880578041077
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.181048601865768
    Y = 0.104528464376926
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.360113620758057
    Y = 0.207911685109138
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.53523313999176
    Y = 0.30901700258255
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.7044885158538819
    Y = 0.406736642122269
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.866025388240814
    Y = 0.5
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.01807391643524
    Y = 0.587785243988037
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.15896821022034
    Y = 0.669130623340607
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.28716456890106
    Y = 0.743144810199738
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.40125858783722
    Y = 0.80901700258255
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.5
    Y = 0.866025388240814
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58230710029602
    Y = 0.913545429706573
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.64727818965912
    Y = 0.95105654001236
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6942013502121
    Y = 0.978147625923157
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.72256243228912
    Y = 0.994521915912628
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73205077648163
    Y = 1
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.72256243228912
    Y = 0.994521915912628
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6942013502121
    Y = 0.978147625923157
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.64727818965912
    Y = 0.95105654001236
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58230710029602
    Y = 0.913545429706573
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.5
    Y = 0.866025388240814
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.40125858783722
    Y = 0.80901700258255
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.28716456890106
    Y = 0.743144810199738
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.15896821022034
    Y = 0.669130623340607
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.01807391643524
    Y = 0.587785243988037
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.866025388240814
    Y = 0.5
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.7044885158538819
    Y = 0.406736642122269
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.53523313999176
    Y = 0.30901700258255
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.360113620758057
    Y = 0.207911685109138
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.181048601865768
    Y = 0.104528464376926
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -2.12115049442606E-16
    Y = 1.22464685258517E-16
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.181048601865768
    Y = -0.104528464376926
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.360113620758057
    Y = -0.207911685109138
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.53523313999176
    Y = -0.30901700258255
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.7044885158538819
    Y = -0.406736642122269
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.866025388240814
    Y = -0.5
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.01807391643524
    Y = -0.587785243988037
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.15896821022034
    Y = -0.669130623340607
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.28716456890106
    Y = -0.743144810199738
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.40125858783722
    Y = -0.80901700258255
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.5
    Y = -0.866025388240814
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58230710029602
    Y = -0.913545429706573
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.64727818965912
    Y = -0.95105654001236
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6942013502121
    Y = -0.978147625923157
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.72256243228912
    Y = -0.994521915912628
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73205077648163
    Y = -1
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.72256243228912
    Y = -0.994521915912628
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6942013502121
    Y = -0.978147625923157
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.64727818965912
    Y = -0.95105654001236
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58230710029602
    Y = -0.913545429706573
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.5
    Y = -0.866025388240814
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.40125858783722
    Y = -0.80901700258255
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.28716456890106
    Y = -0.743144810199738
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.15896821022034
    Y = -0.669130623340607
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.01807391643524
    Y = -0.587785243988037
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.866025388240814
    Y = -0.5
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.7044885158538819
    Y = -0.406736642122269
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.53523313999176
    Y = -0.30901700258255
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.360113620758057
    Y = -0.207911685109138
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.181048601865768
    Y = -0.104528464376926
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = 0.0850311145186424
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.37987357378006
    Y = 0.169130608439445
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.564602136611938
    Y = 0.251377075910568
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = 0.330869406461716
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.913545429706573
    Y = 0.406736642122269
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.07393705844879
    Y = 0.478147596120834
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.22256243228912
    Y = 0.544319868087769
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.35779321193695
    Y = 0.604528486728668
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = 0.658113718032837
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58230710029602
    Y = 0.7044885158538819
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.66913056373596
    Y = 0.743144810199738
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73766672611237
    Y = 0.77365905046463
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.78716456890106
    Y = 0.795696914196014
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.81708192825317
    Y = 0.80901700258255
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.82709085941315
    Y = 0.813473284244537
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.81708192825317
    Y = 0.80901700258255
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.78716456890106
    Y = 0.795696914196014
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73766672611237
    Y = 0.77365905046463
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.66913056373596
    Y = 0.743144810199738
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58230710029602
    Y = 0.7044885158538819
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = 0.658113718032837
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.35779321193695
    Y = 0.604528486728668
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.22256243228912
    Y = 0.544319868087769
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.07393705844879
    Y = 0.478147596120834
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.913545429706573
    Y = 0.406736642122269
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.743144810199738
    Y = 0.330869406461716
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.564602136611938
    Y = 0.251377075910568
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.37987357378006
    Y = 0.169130608439445
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.190983012318611
    Y = 0.0850311145186424
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -2.23754102701048E-16
    Y = 9.96217485357301E-17
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = -0.0850311145186424
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.37987357378006
    Y = -0.169130608439445
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.564602136611938
    Y = -0.251377075910568
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = -0.330869406461716
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.913545429706573
    Y = -0.406736642122269
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.07393705844879
    Y = -0.478147596120834
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.22256243228912
    Y = -0.544319868087769
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.35779321193695
    Y = -0.604528486728668
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = -0.658113718032837
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58230710029602
    Y = -0.7044885158538819
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.66913056373596
    Y = -0.743144810199738
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73766672611237
    Y = -0.77365905046463
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.78716456890106
    Y = -0.795696914196014
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.81708192825317
    Y = -0.80901700258255
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.82709085941315
    Y = -0.813473284244537
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.81708192825317
    Y = -0.80901700258255
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.78716456890106
    Y = -0.795696914196014
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73766672611237
    Y = -0.77365905046463
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.66913056373596
    Y = -0.743144810199738
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58230710029602
    Y = -0.7044885158538819
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = -0.658113718032837
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.35779321193695
    Y = -0.604528486728668
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.22256243228912
    Y = -0.544319868087769
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.07393705844879
    Y = -0.478147596120834
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.913545429706573
    Y = -0.406736642122269
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.743144810199738
    Y = -0.330869406461716
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.564602136611938
    Y = -0.251377075910568
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.37987357378006
    Y = -0.169130608439445
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.190983012318611
    Y = -0.0850311145186424
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.19882495701313
    Y = 0.0646021440625191
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.395471543073654
    Y = 0.128496497869492
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = 0.190983012318611
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.77365905046463
    Y = 0.251377075910568
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = 0.30901700258255
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.1180340051651
    Y = 0.363271266222
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.27276206016541
    Y = 0.413545459508896
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.41354548931122
    Y = 0.459288746118546
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.53884172439575
    Y = 0.5
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.64727818965912
    Y = 0.53523313999176
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73766672611237
    Y = 0.564602136611938
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.80901694297791
    Y = 0.587785243988037
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.86054730415344
    Y = 0.604528486728668
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.89169311523438
    Y = 0.614648342132568
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.90211308002472
    Y = 0.6180340051651
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.89169311523438
    Y = 0.614648342132568
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.86054730415344
    Y = 0.604528486728668
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.80901694297791
    Y = 0.587785243988037
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.73766672611237
    Y = 0.564602136611938
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.64727818965912
    Y = 0.53523313999176
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.53884172439575
    Y = 0.5
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.41354548931122
    Y = 0.459288746118546
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.27276206016541
    Y = 0.413545459508896
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.1180340051651
    Y = 0.363271266222
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.95105654001236
    Y = 0.30901700258255
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.77365905046463
    Y = 0.251377075910568
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.587785243988037
    Y = 0.190983012318611
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.395471543073654
    Y = 0.128496497869492
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.19882495701313
    Y = 0.0646021440625191
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -2.32941657321678E-16
    Y = 7.56873370685124E-17
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.19882495701313
    Y = -0.0646021440625191
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.395471543073654
    Y = -0.128496497869492
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = -0.190983012318611
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.77365905046463
    Y = -0.251377075910568
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = -0.30901700258255
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.1180340051651
    Y = -0.363271266222
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.27276206016541
    Y = -0.413545459508896
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.41354548931122
    Y = -0.459288746118546
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.53884172439575
    Y = -0.5
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.64727818965912
    Y = -0.53523313999176
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73766672611237
    Y = -0.564602136611938
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.80901694297791
    Y = -0.587785243988037
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.86054730415344
    Y = -0.604528486728668
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.89169311523438
    Y = -0.614648342132568
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.90211308002472
    Y = -0.6180340051651
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.89169311523438
    Y = -0.614648342132568
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.86054730415344
    Y = -0.604528486728668
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.80901694297791
    Y = -0.587785243988037
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.73766672611237
    Y = -0.564602136611938
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.64727818965912
    Y = -0.53523313999176
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.53884172439575
    Y = -0.5
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.41354548931122
    Y = -0.459288746118546
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.27276206016541
    Y = -0.413545459508896
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.1180340051651
    Y = -0.363271266222
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.95105654001236
    Y = -0.30901700258255
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.77365905046463
    Y = -0.251377075910568
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.587785243988037
    Y = -0.190983012318611
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.395471543073654
    Y = -0.128496497869492
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.19882495701313
    Y = -0.0646021440625191
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.204488530755043
    Y = 0.0434653796255589
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = 0.0864545404911041
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = 0.128496497869492
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.795696914196014
    Y = 0.169130608439445
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.978147625923157
    Y = 0.207911685109138
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.14988148212433
    Y = 0.24441485106945
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = 0.278240144252777
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.4538106918335
    Y = 0.30901700258255
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58267605304718
    Y = 0.336408168077469
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6942013502121
    Y = 0.360113620758057
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.78716456890106
    Y = 0.37987357378006
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.86054730415344
    Y = 0.395471543073654
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.91354548931122
    Y = 0.406736642122269
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.94557845592499
    Y = 0.413545459508896
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.95629525184631
    Y = 0.415823370218277
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.94557845592499
    Y = 0.413545459508896
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.91354548931122
    Y = 0.406736642122269
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.86054730415344
    Y = 0.395471543073654
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.78716456890106
    Y = 0.37987357378006
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.6942013502121
    Y = 0.360113620758057
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.58267605304718
    Y = 0.336408168077469
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.4538106918335
    Y = 0.30901700258255
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.30901694297791
    Y = 0.278240144252777
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.14988148212433
    Y = 0.24441485106945
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.978147625923157
    Y = 0.207911685109138
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.795696914196014
    Y = 0.169130608439445
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.604528486728668
    Y = 0.128496497869492
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.406736642122269
    Y = 0.0864545404911041
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.204488530755043
    Y = 0.0434653796255589
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -2.39577075202452E-16
    Y = 5.09236770372102E-17
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.204488530755043
    Y = -0.0434653796255589
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = -0.0864545404911041
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = -0.128496497869492
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.795696914196014
    Y = -0.169130608439445
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.978147625923157
    Y = -0.207911685109138
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.14988148212433
    Y = -0.24441485106945
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = -0.278240144252777
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.4538106918335
    Y = -0.30901700258255
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58267605304718
    Y = -0.336408168077469
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6942013502121
    Y = -0.360113620758057
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.78716456890106
    Y = -0.37987357378006
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.86054730415344
    Y = -0.395471543073654
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.91354548931122
    Y = -0.406736642122269
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.94557845592499
    Y = -0.413545459508896
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.95629525184631
    Y = -0.415823370218277
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.94557845592499
    Y = -0.413545459508896
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.91354548931122
    Y = -0.406736642122269
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.86054730415344
    Y = -0.395471543073654
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.78716456890106
    Y = -0.37987357378006
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.6942013502121
    Y = -0.360113620758057
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.58267605304718
    Y = -0.336408168077469
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.4538106918335
    Y = -0.30901700258255
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.30901694297791
    Y = -0.278240144252777
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.14988148212433
    Y = -0.24441485106945
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.978147625923157
    Y = -0.207911685109138
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.795696914196014
    Y = -0.169130608439445
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.604528486728668
    Y = -0.128496497869492
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.406736642122269
    Y = -0.0864545404911041
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.204488530755043
    Y = -0.0434653796255589
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    Z = 2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = 0.0218524001538754
    Z = 1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = 0.0434653796255589
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.614648342132568
    Y = 0.0646021440625191
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = 0.0850311145186424
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = 0.104528464376926
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.16913056373596
    Y = 0.122880578041077
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33093011379242
    Y = 0.139886394143105
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = 0.1553595662117
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.60917019844055
    Y = 0.169130608439445
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.72256243228912
    Y = 0.181048601865768
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.81708192825317
    Y = 0.190983012318611
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.89169311523438
    Y = 0.19882495701313
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.94557845592499
    Y = 0.204488530755043
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.97814762592316
    Y = 0.207911685109138
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.98904383182526
    Y = 0.209056928753853
    Z = 1.22464685258517E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.97814762592316
    Y = 0.207911685109138
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.94557845592499
    Y = 0.204488530755043
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.89169311523438
    Y = 0.19882495701313
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.81708192825317
    Y = 0.190983012318611
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.72256243228912
    Y = 0.181048601865768
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.60917019844055
    Y = 0.169130608439445
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.47814762592316
    Y = 0.1553595662117
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.33093011379242
    Y = 0.139886394143105
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -1.16913056373596
    Y = 0.122880578041077
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.994521915912628
    Y = 0.104528464376926
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.80901700258255
    Y = 0.0850311145186424
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.614648342132568
    Y = 0.0646021440625191
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.413545459508896
    Y = 0.0434653796255589
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -0.207911685109138
    Y = 0.0218524001538754
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = -2.43587617389024E-16
    Y = 2.56020903001859E-17
    Z = -2
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = -0.0218524001538754
    Z = -1.98904383182526
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = -0.0434653796255589
    Z = -1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.614648342132568
    Y = -0.0646021440625191
    Z = -1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = -0.0850311145186424
    Z = -1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = -0.104528464376926
    Z = -1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.16913056373596
    Y = -0.122880578041077
    Z = -1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33093011379242
    Y = -0.139886394143105
    Z = -1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = -0.1553595662117
    Z = -1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.60917019844055
    Y = -0.169130608439445
    Z = -1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.72256243228912
    Y = -0.181048601865768
    Z = -1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.81708192825317
    Y = -0.190983012318611
    Z = -0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.89169311523438
    Y = -0.19882495701313
    Z = -0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.94557845592499
    Y = -0.204488530755043
    Z = -0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.97814762592316
    Y = -0.207911685109138
    Z = -0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.98904383182526
    Y = -0.209056928753853
    Z = -3.67394029305771E-16
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.97814762592316
    Y = -0.207911685109138
    Z = 0.209056928753853
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.94557845592499
    Y = -0.204488530755043
    Z = 0.415823370218277
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.89169311523438
    Y = -0.19882495701313
    Z = 0.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.81708192825317
    Y = -0.190983012318611
    Z = 0.813473284244537
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.72256243228912
    Y = -0.181048601865768
    Z = 1
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.60917019844055
    Y = -0.169130608439445
    Z = 1.17557048797607
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.47814762592316
    Y = -0.1553595662117
    Z = 1.33826124668121
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.33093011379242
    Y = -0.139886394143105
    Z = 1.48628962039948
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 1.16913056373596
    Y = -0.122880578041077
    Z = 1.6180340051651
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.994521915912628
    Y = -0.104528464376926
    Z = 1.73205077648163
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.80901700258255
    Y = -0.0850311145186424
    Z = 1.82709085941315
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.614648342132568
    Y = -0.0646021440625191
    Z = 1.90211308002472
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.413545459508896
    Y = -0.0434653796255589
    Z = 1.95629525184631
    }
$!Streamtrace Add
  StreamType = VolumeLine
  StartPos
    {
    X = 0.207911685109138
    Y = -0.0218524001538754
    Z = 1.98904383182526
    }
$!FrameSetup Default3DStreamtrace{StreamType = VolumeLine}
$!StreamAttributes CellFraction = 1
$!StreamAttributes MaxSteps = 750
$!GlobalThreeDVector UVar = 5
$!GlobalThreeDVector VVar = 6
$!GlobalThreeDVector WVar = 7
$!StreamtraceLayers Show = No
$!StreamtraceLayers Show = Yes
$!ExtractStreamtraces 
$!WriteDataSet  "|STREAMFILE|"
  IncludeText = No
  IncludeGeom = No
  IncludeCustomLabels = No
  IncludeDataShareLinkage = Yes
  ZoneList =  [2-|NUMZONES|]
  VarPositionList =  [1-3,5-7,35-38,44-52]
  Binary = No
  UsePointFormat = Yes
  Precision = 9
  TecplotVersionToWrite = TecplotCurrent
