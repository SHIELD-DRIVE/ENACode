#!MC 1410
# Step 2b: interpolates the neutral variables (15-39) from zone 2 onto the grid (zone 1), then
# deletes zone 2.
$!LinearInterpolate 
  SourceZones =  [2]
  DestinationZone = 1
  VarList =  [15-39]
  LinearInterPConst = 0
  LinearInterpMode = DontChange
$!DeleteZones  [2]
