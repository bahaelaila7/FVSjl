data/fire/photo_series.csv is produced from the oracle source, not transcribed:
  B=/workspace/ForestVegetationSimulator/bin/FVSne_buildDir   (fmphotoval.f is byte-identical in the CS/LS/NE/SN/ON/west builds)
  cp $B/fmphotoval.f $B/*.F77 . && gfortran-16 -std=legacy -w -O0 -fno-automatic -o genphoto gen.f fmphotoval.f stubs.f && ./genphoto > photo_ne.csv
(gen.f calls FMPHOTOVAL for every reference 1..32 / code 1..200 and prints the valid rows; stubs.f replaces DBCHK/RCDSET.)
