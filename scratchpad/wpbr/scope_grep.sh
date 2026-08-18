#!/bin/bash
FVS=/workspace/ForestVegetationSimulator
echo "===== exbrus.f stub (entry points to provide) ====="
grep -iE "entry|subroutine" "$FVS"/base/exbrus.f
echo; echo "===== keywds.f: BLISTER / wpbr option number ====="
grep -niE "blister|wpbr|brus" "$FVS"/base/keywds.f
echo; echo "===== driver call sites in base (grincr/gradd/fvs/etc) ====="
grep -rniE "call br(in|setp|dam|gi|init|updt|star|atv)\b|blister" "$FVS"/base/*.f | head -40
echo; echo "===== RNG in wpbr (brann.f) ====="
grep -niE "16807|2147483647|dmod|seed|entry" "$FVS"/wpbr/brann.f
echo; echo "===== which variant source lists reference wpbr ====="
grep -rilE "wpbr|brcank|brinit|exbrus" "$FVS"/*/CMakeLists.txt 2>/dev/null | head
grep -rilE "wpbr" "$FVS"/bin/*_buildDir/*.txt 2>/dev/null | head
