# Extract CI large-tree DDS coefficient arrays from ci/dgf.f → data/centralidaho/dg_coeffs_1d.csv.
# 1-D per-species DDS arrays (len 19). 2-D arrays (DGFOR/DGDS/DGCCFA/DGHAB/ICHBCL/MAPLOC/MAP*) are
# extracted separately (per-column DATA forms). Uses tools/fortran_data_extract.grab (repeat-expands n*v).
import sys; sys.path.insert(0, 'tools')
from fortran_data_extract import grab
SRC = "/workspace/ForestVegetationSimulator/bin/FVSci_buildDir/dgf.f"
ONED = ["DGLD","DGCR","DGCRSQ","DGBAL","DGDBAL","DGBA","DGLBA","DGPCCF","DGEL","DGEL2","DGSLOP","DGSLSQ","DGCASP","DGSASP"]
cols = {}
for n in ONED:
    v = grab(SRC, n)
    assert v and len(v) == 19, f"{n}: got {len(v) if v else None}"
    cols[n] = v
with open("data/centralidaho/dg_coeffs_1d.csv", "w") as f:
    f.write("species_index," + ",".join(ONED) + "\n")
    for i in range(19):
        f.write(str(i+1) + "," + ",".join(repr(cols[n][i]) for n in ONED) + "\n")
print("wrote data/centralidaho/dg_coeffs_1d.csv (19 rows x %d DDS arrays)" % len(ONED))
