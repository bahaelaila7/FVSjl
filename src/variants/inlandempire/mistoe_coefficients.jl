# IE MISTOE (dwarf mistletoe) effect coefficients — transcribed from mistoe/misintie.f DATA
# (the "Interior Empire" NI/KT/EM init). These are the EFFECT tables only; the infection (misinf.f) +
# spread (mistoe.f) subsystem is a separate port. IE has NO height-growth reduction from DM (AHGP all 1.0),
# so IE MISTOE affects only diameter growth (IE_MIS_DGP) and mortality (IE_MIS_PMC).
#
# DMR = dwarf-mistletoe rating 0..6 per tree (IMIST). Affected species (MISFIT): sp2 WL(larch), 3 DF,
# 7 LP, 10 PP, 12 WB, 13 LM. WB/LM reuse LP's curves; WL reuses DF's.

# MISFIT — 1 if species is DM-affected (mistoe/misintie.f AFIT)
const IE_MIS_FIT = Int32[0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# DGPDMR — diameter-growth MULTIPLIER by species (col) × DMR 0..6 (row). Access IE_MIS_DGP[DMR+1, sp].
# (misintie.f ADGP(MAXSP,7); flat is per-species DMR-fastest ⇒ reshape 7×23 column-major.)
const IE_MIS_DGP = reshape(Float32[
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp1
    1.0,0.94,0.92,0.88,0.84,0.58,0.54,    # sp2 WL/larch
    1.0,0.98,0.97,0.85,0.80,0.52,0.44,    # sp3 DF
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp4
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp5
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp6
    1.0,1.0,1.0,1.0,0.94,0.80,0.59,       # sp7 LP
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp8
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp9
    1.0,1.0,1.0,0.98,0.86,0.73,0.50,      # sp10 PP
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp11
    1.0,1.0,1.0,1.0,0.94,0.80,0.59,       # sp12 WB (=LP)
    1.0,1.0,1.0,1.0,0.94,0.80,0.59,       # sp13 LM (=LP)
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp14
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp15
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp16
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp17
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp18
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp19
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp20
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp21
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp22
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,          # sp23
], 7, 23)

# PMCSP — DM-mortality coefficients, 3 per species (misintie.f APMC). Access IE_MIS_PMC[1:3, sp].
const IE_MIS_PMC = reshape(Float32[
    0.0,0.0,0.0,                          # sp1
    0.01319,-0.01627,0.00822,             # sp2 WL (=DF)
    0.01319,-0.01627,0.00822,             # sp3 DF
    0.0,0.0,0.0,                          # sp4
    0.0,0.0,0.0,                          # sp5
    0.0,0.0,0.0,                          # sp6
    0.00112,0.02170,-0.00171,             # sp7 LP
    0.0,0.0,0.0,                          # sp8
    0.0,0.0,0.0,                          # sp9
    0.00681,-0.00580,0.00935,             # sp10 PP
    0.0,0.0,0.0,                          # sp11
    0.00112,0.02170,-0.00171,             # sp12 WB (=LP)
    0.00112,0.02170,-0.00171,             # sp13 LM (=LP)
    0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0,   # sp14-18
    0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0,   # sp19-23
], 3, 23)
