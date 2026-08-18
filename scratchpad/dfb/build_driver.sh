#!/bin/bash
set -e
D=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb
BD=/workspace/ForestVegetationSimulator/bin/FVSbc_buildDir
cd "$D"
gfortran-16 -std=legacy -w -fno-automatic -I"$BD" -I/workspace/ForestVegetationSimulator/dfb/common \
  driver_dfb.f obj/dfber.o obj/dfbprb.o obj/dfbdbh.o obj/dfbind.o obj/dfblkdbc.o -o driver_dfb
echo "built driver_dfb"
