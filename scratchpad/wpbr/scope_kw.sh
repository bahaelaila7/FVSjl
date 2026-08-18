#!/bin/bash
FVS=/workspace/ForestVegetationSimulator
echo "===== keywds.f array context + BRUST dispatch ====="
grep -n "BRUST" "$FVS"/base/keywds.f
echo "--- how keywds maps name->number (GOTO / KOItabl) ---"
sed -n '20,120p' "$FVS"/base/keywds.f
