#!MC 1410
# Step 2a of the single-ion input preparation (tools/tecplot_prep.py).
# Starts from: the spherical ENA grid, loaded by the driver (input "grid"; zone 1).
# Does: appends the fixed neutrals |NEUTRALFIXED| (.szplt, from step 1) as zone 2, matched by
#       variable name.
# Steps 2a-2d are interpolating_files.mcr split in four, so that each part is timed on its own;
# together they do exactly what the recorded macro did.
# Read with the subzone loader, as recorded in Tecplot 2025 R2 (interpolating_files_2025.mcr).
$!ReadDataSet  '"STANDARDSYNTAX" "1.0" "FILENAME_FILE" "|NEUTRALFIXED|"'
  DataSetReader = 'Tecplot Subzone Data Loader'
  VarNameList = '"X AU";"X" "Y AU";"Y" "Z AU";"Z" "Rho amu/cm3" "U_x km/s" "U_y km/s" "U_z km/s" "B_x nT" "B_y nT" "B_z nT" "P dyne/cm^2" "J_x uA/m^2" "J_y uA/m^2" "J_z uA/m^2" "nH_SSW" "UxH_SSW";"Ux_SSW" "UyH_SSW";"Uy_SSW" "UzH_SSW";"Uz_SSW" "TH_SSW";"T_SSW" "nH_IHS" "UxH_IHS";"Ux_IHS" "UyH_IHS";"Uy_IHS" "UzH_IHS";"Uz_IHS" "TH_IHS";"T_IHS" "nH_OHS" "UxH_OHS";"Ux_OHS" "UyH_OHS";"Uy_OHS" "UzH_OHS";"Uz_OHS" "TH_OHS";"T_OHS" "nH_ISM" "UxH_ISM";"Ux_ISM" "UyH_ISM";"Uy_ISM" "UzH_ISM";"Uz_ISM" "TH_ISM";"T_ISM" "nH" "UxH" "UyH" "UzH" "TH" "Hplim amu/cm3"'
  ReadDataOption = Append
  ResetStyle = No
  AssignStrandIDs = No
  InitialPlotType = Automatic
  InitialPlotFirstZoneOnly = No
  AddZonesToExistingStrands = No
  VarLoadMode = ByName
