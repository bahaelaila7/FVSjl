#!/bin/bash
set -e
BD=/workspace/ForestVegetationSimulator/bin/FVSon_buildDir
OW=/workspace/.onwork
FLAGS="-std=legacy -w -fno-automatic -fPIC -c -I$BD"
gfortran-16 $FLAGS -o $OW/dbs_fiavbc_trls.o   /workspace/ForestVegetationSimulator/dbsqlite/dbs_fiavbc_trls.f
gfortran-16 $FLAGS -o $OW/dbs_fiavbc_cutlst.o /workspace/ForestVegetationSimulator/dbsqlite/dbs_fiavbc_cutlst.f
gfortran-16 $FLAGS -o $OW/dbs_fiavbc_atrtls.o /workspace/ForestVegetationSimulator/dbsqlite/dbs_fiavbc_atrtls.f
gfortran-16 $FLAGS -o $OW/dbsreference.o      /workspace/ForestVegetationSimulator/vdbsqlite/dbsreference.f
echo "compiled 4 missing objects"
