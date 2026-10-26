#!MC 1410
# Step 2d: interpolates the plasma variables (4-14 and 40) from zone 2 onto the grid (zone 1),
# then deletes zone 2.
$!LinearInterpolate 
  SourceZones =  [2]
  DestinationZone = 1
  VarList =  [4-14,40]
  LinearInterPConst = 0
  LinearInterpMode = DontChange
$!DeleteZones  [2]
