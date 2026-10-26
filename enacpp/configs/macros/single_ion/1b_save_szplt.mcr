#!MC 1410
# Step 1b of the single-ion input preparation (configs/prep_single_ion.toml).
# Starts from: the fixed neutrals written by ena_fixneutrals (|NEUTRALPLT|, a binary .plt),
#       loaded by the driver (work file "neutral_plt").
# Does: nothing to the data; saves it as .szplt, with the command Tecplot 2025 R2 records for
#       File > Write Data File > .szplt. Tecplot 2025 R2 interpolates quickly (step 2b) only from a
#       .szplt neutral file.
# Writes: |NEUTRALFIXED| (an intermediate .szplt file in the work folder), read by step 2a.
$!ExtendedCommand 
  CommandProcessorID = 'Tecplot Subzone Data Tools'
  Command = 'WRITEDATASET FILENAME="|NEUTRALFIXED|"'
