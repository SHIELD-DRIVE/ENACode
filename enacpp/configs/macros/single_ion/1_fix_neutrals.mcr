#!MC 1410
# Step 1 of the single-ion input preparation in Tecplot (configs/prep_single_ion_tecplot.toml).
# The standard route (prep_single_ion.toml) does the same with the C++ program ena_fixneutrals.
# Starts from: the neutral solution of this snapshot, loaded by the driver (input "neutral_source").
# Does: renames the population variables, computes the total nH, UxH, UyH, UzH, TH and the
#       temperature of each population.
# Writes: |NEUTRALFIXED| (an intermediate .szplt file in the work folder), read by step 2a.
# Generalised from fixing_neutrals.mcr: the output file name was changed, and the file is
# written as .szplt instead of ASCII.
$!RenameDataSetVar 
  Var = 1
  Name = 'X AU'
$!RenameDataSetVar 
  Var = 2
  Name = 'Y AU'
$!RenameDataSetVar 
  Var = 3
  Name = 'Z AU'
$!RenameDataSetVar 
  Var = 10
  Name = 'nH_OHS'
$!RenameDataSetVar 
  Var = 12
  Name = 'UyH_OHS'
$!RenameDataSetVar 
  Var = 11
  Name = 'UxH_OHS'
$!RenameDataSetVar 
  Var = 13
  Name = 'UzH_OHS'
$!RenameDataSetVar 
  Var = 21
  Name = 'nH_IHS'
$!RenameDataSetVar 
  Var = 22
  Name = 'UxH_IHS'
$!RenameDataSetVar 
  Var = 23
  Name = 'UyH_IHS'
$!RenameDataSetVar 
  Var = 24
  Name = 'UzH_IHS'
$!RenameDataSetVar 
  Var = 32
  Name = 'nH_SSW'
$!RenameDataSetVar 
  Var = 33
  Name = 'UxH_SSW'
$!RenameDataSetVar 
  Var = 34
  Name = 'UyH_SSW'
$!RenameDataSetVar 
  Var = 35
  Name = 'UzH_SSW'
$!RenameDataSetVar 
  Var = 43
  Name = 'nH_ISM'
$!RenameDataSetVar 
  Var = 44
  Name = 'UxH_ISM'
$!RenameDataSetVar 
  Var = 45
  Name = 'UyH_ISM'
$!RenameDataSetVar 
  Var = 46
  Name = 'UzH_ISM'
$!AlterData 
  Equation = '{nH}=v10+v21+v32+v43'
$!AlterData 
  Equation = '{UxH}=(v10*v11+v21*v22+v32*v33+v43*v44)/{nH}'
$!AlterData 
  Equation = '{UyH}=(v10*v12+v21*v23+v32*v34+v43*v45)/{nH}'
$!AlterData 
  Equation = '{UzH}=(v10*v13+v21*v24+v32*v35+v43*v46)/{nH}'
$!AlterData 
  Equation = '{Pxx}=(v14+v25+v36+v47)*1e-9+(v10*(v11-{UxH})**2+v21*(v22-{UxH})**2+v32*(v33-{UxH})**2+v43*(v44-{UxH})**2)*1.67e-27*(1000**2)*(100**3)'
$!AlterData 
  Equation = '{Pyy}=(v15+v26+v37+v48)*1e-9+(v10*(v12-{UyH})**2+v21*(v23-{UyH})**2+v32*(v34-{UyH})**2+v43*(v45-{UyH})**2)*1.67e-27*(1000**2)*(100**3)'
$!AlterData 
  Equation = '{Pzz}=(v16+v27+v38+v49)*1e-9+(v10*(v13-{UzH})**2+v21*(v24-{UzH})**2+v32*(v35-{UzH})**2+v43*(v46-{UzH})**2)*1.67e-27*(1000**2)*(100**3)'
$!AlterData 
  Equation = '{P}=1/3*({Pxx}+{Pyy}+{Pzz})'
$!AlterData 
  Equation = '{TH}={P}/({nH}*(100**3)*1.38e-23+1e-30)'
$!AlterData 
  Equation = '{P_OHS}=1/3*({pxxPop1}+{pyyPop1}+{pzzPop1})*1e-9'
$!AlterData 
  Equation = '{TH_OHS}={P_OHS}/({nH_OHS}*(100**3)*1.38e-23+1e-30)'
$!AlterData 
  Equation = '{P_IHS}=1/3*({pxxPop2}+{pyyPop2}+{pzzPop2})*1e-9'
$!AlterData 
  Equation = '{TH_IHS}={P_IHS}/({nH_IHS}*(100**3)*1.38e-23+1e-30)'
$!AlterData 
  Equation = '{P_SSW}=1/3*({pxxPop3}+{pyyPop3}+{pzzPop3})*1e-9'
$!AlterData 
  Equation = '{TH_SSW}={P_SSW}/({nH_SSW}*(100**3)*1.38e-23+1e-30)'
$!AlterData 
  Equation = '{P_ISM}=1/3*({pxxPop4}+{pyyPop4}+{pzzPop4})*1e-9'
$!AlterData 
  Equation = '{TH_ISM}={P_ISM}/({nH_ISM}*(100**3)*1.38e-23+1e-30)'
$!DeleteVars  [4-9]
$!DeleteVars  [8-14]
$!DeleteVars  [12-18]
$!DeleteVars  [16-22]
$!DeleteVars  [20-26]
$!DeleteVars  [24-26]
# Saved as .szplt (Tecplot's subzone format), with the command Tecplot records for File > Write
# Data File > .szplt. Tecplot 2025 R2 interpolates quickly from a .szplt neutral file in step 2b,
# but stalls for hours with a binary .plt or ASCII .dat one (2023 R2 was as slow with .plt).
$!ExtendedCommand 
  CommandProcessorID = 'Tecplot Subzone Data Tools'
  Command = 'WRITEDATASET FILENAME="|NEUTRALFIXED|"'
