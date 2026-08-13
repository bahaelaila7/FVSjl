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

# ---------------------------------------------------------------------------
# UTAH 24-species mistletoe tables (misintut.f AFIT/ADGP/APMC). UT's species order
# DIFFERS from IE — DF(3)/LP(7)/PP(10) align, but UT's pinyon/woodland hosts PI(11)/
# PM(14)/GB(17) sit at UT-specific indices absent from the IE table. jl previously
# applied the IE table to UT ⇒ UT pinyon (FIA133→sp14) got mortality-rate 0 while live
# kills it (misintut.f AFIT(14)=1, surrogate LP coefficients) ⇒ UT over-grew mistletoe
# stands. LP=(.00112,.02170,-.00171) DF=(.01319,-.01627,.00822) GF/WF=(0,.00159,.00508)
# PP=(.00681,-.00580,.00935). Surrogates per misintut.f comments (PM→LP, GB→PP, …).
const UT_MIS_FIT = Int32[1,1,1,1,1,0,1,1,1,1, 1,0,0,1,0,0,1,0,0,0, 0,0,0,0]
const UT_MIS_PMC = reshape(Float32[
    0.00112,0.02170,-0.00171,   # 1 WB(LP)
    0.00112,0.02170,-0.00171,   # 2 LM(LP)
    0.01319,-0.01627,0.00822,   # 3 DF
    0.0,0.00159,0.00508,        # 4 WF(GF)
    0.00681,-0.00580,0.00935,   # 5 BS(PP)
    0.0,0.0,0.0,                # 6 AS
    0.00112,0.02170,-0.00171,   # 7 LP
    0.01319,-0.01627,0.00822,   # 8 ES(DF)
    0.0,0.00159,0.00508,        # 9 AF(GF)
    0.00681,-0.00580,0.00935,   # 10 PP
    0.00112,0.02170,-0.00171,   # 11 PI(LP)
    0.0,0.0,0.0,                # 12 WJ
    0.0,0.0,0.0,                # 13 GO
    0.00112,0.02170,-0.00171,   # 14 PM(LP)
    0.0,0.0,0.0,                # 15 RM
    0.0,0.0,0.0,                # 16 UJ
    0.00681,-0.00580,0.00935,   # 17 GB(PP)
    0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0,   # 18 NC, 19 FC, 20 MC
    0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0,   # 21 BI, 22 BE, 23 OS, 24 OH
], 3, 24)
const UT_MIS_DGP = reshape(Float32[
    1.0,1.0,1.0,1.0,0.94,0.80,0.59,   # 1 WB(LP)
    1.0,1.0,1.0,1.0,0.94,0.80,0.59,   # 2 LM(LP)
    1.0,0.98,0.97,0.85,0.80,0.52,0.44,# 3 DF
    1.0,1.0,1.0,0.98,0.95,0.70,0.50,  # 4 WF(GF)
    1.0,1.0,1.0,0.98,0.86,0.73,0.50,  # 5 BS(PP)
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,      # 6 AS
    1.0,1.0,1.0,1.0,0.94,0.80,0.59,   # 7 LP
    1.0,0.98,0.97,0.85,0.80,0.52,0.44,# 8 ES(DF)
    1.0,1.0,1.0,0.98,0.95,0.70,0.50,  # 9 AF(GF)
    1.0,1.0,1.0,0.98,0.86,0.73,0.50,  # 10 PP
    1.0,1.0,1.0,1.0,0.94,0.80,0.59,   # 11 PI(LP)
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,      # 12 WJ
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,      # 13 GO
    1.0,1.0,1.0,1.0,0.94,0.80,0.59,   # 14 PM(LP)
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,      # 15 RM
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,      # 16 UJ
    1.0,1.0,1.0,0.98,0.86,0.73,0.50,  # 17 GB(PP)
    1.0,1.0,1.0,1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0,1.0,1.0,1.0,  # 18-20
    1.0,1.0,1.0,1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0,1.0,1.0,1.0,  # 21-24
], 7, 24)

# ---------------------------------------------------------------------------
# EM/KT/BM/TT/CI per-variant mistletoe tables — parser-extracted verbatim from each
# misint{v}.f (AFIT/ADGP/APMC), Julia-normalized floats. Each variant's species ORDER
# differs; DF/LP/PP align with IE (validated) but variant-specific hosts (white pine,
# grand fir, woodland) need their own coef. Extraction verified against IE/UT.
const EM_MIS_FIT = Int32[1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
const EM_MIS_PMC = reshape(Float32[
    0.00112, 0.02170, -0.00171,
    0.0, 0.0, 0.0,
    0.01319, -0.01627, 0.00822,
    0.00112, 0.02170, -0.00171,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.00112, 0.02170, -0.00171,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0
], 3, 19)
const EM_MIS_DGP = reshape(Float32[
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 0.98, 0.97, 0.85, 0.80, 0.52, 0.44,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
], 7, 19)

const KT_MIS_FIT = Int32[0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0]
const KT_MIS_PMC = reshape(Float32[
    0.0, 0.0, 0.0,
    0.01319, -0.01627, 0.00822,
    0.01319, -0.01627, 0.00822,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.00112, 0.02170, -0.00171,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.00681, -0.00580, 0.00935,
    0.0, 0.0, 0.0
], 3, 11)
const KT_MIS_DGP = reshape(Float32[
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 0.94, 0.92, 0.88, 0.84, 0.58, 0.54,
    1.0, 0.98, 0.97, 0.85, 0.80, 0.52, 0.44,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 0.98, 0.86, 0.73, 0.50,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
], 7, 11)

const BM_MIS_FIT = Int32[0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]
const BM_MIS_PMC = reshape(Float32[
    0.0, 0.0, 0.0,
    0.01319, -0.01627, 0.00822,
    0.01319, -0.01627, 0.00822,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.00112, 0.02170, -0.00171,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.00681, -0.00580, 0.00935,
    0.00112, 0.02170, -0.00171,
    0.00112, 0.02170, -0.00171,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0
], 3, 18)
const BM_MIS_DGP = reshape(Float32[
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 0.94, 0.92, 0.88, 0.84, 0.58, 0.54,
    1.0, 0.98, 0.97, 0.85, 0.80, 0.52, 0.44,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 0.98, 0.86, 0.73, 0.50,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
], 7, 18)

const TT_MIS_FIT = Int32[1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
const TT_MIS_PMC = reshape(Float32[
    0.00112, 0.02170, -0.00171,
    0.00112, 0.02170, -0.00171,
    0.01319, -0.01627, 0.00822,
    0.00112, 0.02170, -0.00171,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.00112, 0.02170, -0.00171,
    0.0, 0.0, 0.0,
    0.0, 0.00159, 0.00508,
    0.00681, -0.00580, 0.00935,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0
], 3, 18)
const TT_MIS_DGP = reshape(Float32[
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 0.98, 0.97, 0.85, 0.80, 0.52, 0.44,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 0.98, 0.95, 0.70, 0.50,
    1.0, 1.0, 1.0, 0.98, 0.86, 0.73, 0.50,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
], 7, 18)

const CI_MIS_FIT = Int32[1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0]
const CI_MIS_PMC = reshape(Float32[
    0.00112, 0.02170, -0.00171,
    0.01319, -0.01627, 0.00822,
    0.01319, -0.01627, 0.00822,
    0.0, 0.00159, 0.00508,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.00112, 0.02170, -0.00171,
    0.0, 0.0, 0.0,
    0.0, 0.00159, 0.00508,
    0.00681, -0.00580, 0.00935,
    0.00112, 0.02170, -0.00171,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.00112, 0.02170, -0.00171,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0
], 3, 19)
const CI_MIS_DGP = reshape(Float32[
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 0.94, 0.92, 0.88, 0.84, 0.58, 0.54,
    1.0, 0.98, 0.97, 0.85, 0.80, 0.52, 0.44,
    1.0, 1.0, 1.0, 0.98, 0.95, 0.70, 0.50,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 0.98, 0.95, 0.70, 0.50,
    1.0, 1.0, 1.0, 0.98, 0.86, 0.73, 0.50,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.94, 0.80, 0.59,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
], 7, 19)

# Per-variant mistletoe table dispatch (misint{v}.f DATA differs by species order). Every western
# Wykoff variant now uses its OWN table; IE is the fallback (native). Returns (FIT, DGP, PMC, MAXSP).
@inline function _mis_tables(v)
    v isa Utah          && return (UT_MIS_FIT, UT_MIS_DGP, UT_MIS_PMC, 24)
    v isa EasternMontana&& return (EM_MIS_FIT, EM_MIS_DGP, EM_MIS_PMC, 19)
    v isa Kootenai      && return (KT_MIS_FIT, KT_MIS_DGP, KT_MIS_PMC, 11)
    v isa BlueMountains && return (BM_MIS_FIT, BM_MIS_DGP, BM_MIS_PMC, 18)
    v isa Teton         && return (TT_MIS_FIT, TT_MIS_DGP, TT_MIS_PMC, 18)
    v isa CentralIdaho  && return (CI_MIS_FIT, CI_MIS_DGP, CI_MIS_PMC, 19)
    v isa BritishColumbia && return (BC_MIS_FIT, BC_MIS_DGP, BC_MIS_PMC, 15)  # NEWSPRED C6 payoff (misintbc.f)
    return (IE_MIS_FIT, IE_MIS_DGP, IE_MIS_PMC, 23)   # InlandEmpire (native table)
end

# The DM growth/mortality EFFECTS gate. BC is included here (its DMR comes from NEWSPRED/dm_tregro!,
# not the base ie_mistoe! spread) so the base misdgf/mismrt effects apply to BC's per-tree t.dmr —
# but BC is NOT in _ie_mis_variant, so it does NOT run the base (non-spatial) ie_mistoe! spread.
@inline _dm_effects_variant(v)::Bool = _ie_mis_variant(v) || v isa BritishColumbia

# ============================================================================
# IE MISTOE effect kernels + apply steps (mistoe/misdgf.f + mismrt.f — the SHARED
# equations, IE coefficients). Mirrors the validated CentralRockies path
# (data/centralrockies/dwarf_mistletoe.jl + dwarf_mistletoe_model.jl). DMR is seeded
# generically from tree damage codes 30-34 in treeinput.jl (no IE-specific input needed).
# IE has NO height-growth reduction (AHGP all 1.0) and the INFECTION/SPREAD subsystem
# (misinf/mistoe.f DMR intensification) is a separate later layer — the damage-code path
# gives a STATIC per-tree DMR, sufficient to drive effects + mortality.
# ============================================================================

"misdgf.f: DG diameter-growth multiplier for (species, DMR). 1.0 for DMR 0 / unaffected species."
# The western dwarf-mistletoe EQUATIONS (mistoe.f/misdgf.f/mismrt.f) are shared, but the misint{v}.f species
# DATA (AFIT/ADGP/APMC) is PER-VARIANT — each variant's species ORDER differs. IE's table is correct for a
# variant only where its species index matches IE's (DF=3/LP=7/PP=10 align across the Wykoff cluster). UT's
# pinyon/woodland hosts (PI/PM/GB) sit at UT-specific indices ⇒ UT needs its own table (ported above);
# EM/KT/BM/TT/CI keep IE's for now (correct for their DF/LP/PP; variant woodland hosts = documented follow-up).
# CR has its own 38-sp cr_mistoe!. Gated per-tree DMR ⇒ INERT on stands without dwarf-mistletoe ratings.
@inline _ie_mis_variant(v)::Bool = v isa InlandEmpire || v isa Kootenai || v isa EasternMontana ||
    v isa BlueMountains || v isa Utah || v isa Teton || v isa CentralIdaho

@inline ie_dm_dg_mult(dgp, maxsp::Integer, sp::Integer, dmr::Integer) =
    (sp < 1 || sp > maxsp) ? 1f0 : @inbounds dgp[dmr + 1, sp]

"""
    ie_dm_mortality_rate(pmc, maxsp, sp, dmr, dbh, fint; dmmmlt=1.0) -> Float32

Periodic DM-induced mortality proportion (mismrt.f:155-183): quadratic in DMR, species
multiplier, +20% for DBH<9, floored at 0, capped 0.71(<9)/0.5(≥9), then annualized to the
cycle length. Returns 0 for DMR 0. `pmc`/`maxsp` are the per-variant table (see _mis_tables).
"""
function ie_dm_mortality_rate(pmc, maxsp::Integer, sp::Integer, dmr::Integer, dbh::Real, fint::Real; dmmmlt::Real = 1.0)
    dmr == 0 && return 0.0f0
    (sp < 1 || sp > maxsp) && return 0.0
    b0 = pmc[1, sp]; b1 = pmc[2, sp]; b2 = pmc[3, sp]
    m = b0 + b1 * dmr + b2 * dmr * dmr
    m *= dmmmlt
    small = dbh < 9.0
    small && (m *= 1.2)
    m < 0.0 && (m = 0.0)
    cap = small ? 0.71 : 0.5
    m > cap && (m = cap)
    return Float32(1.0 - (1.0 - m)^(fint / 10.0))
end

"""
    ie_dm_growth_loss!(s, stash)

misdgf.f (applied at dgdriv.f:230, post-DG-driver): multiply each infected tree's diameter
growth (central + tripled dgU/dgL) by IE_MIS_DGP[DMR+1,sp], using START-of-cycle DMR.
No-op for non-IE / uninfected. Deterministic.
"""
function ie_dm_growth_loss!(s::StandState, stash)
    _dm_effects_variant(s.variant) || return
    t = s.trees
    _, dgp, _, maxsp = _mis_tables(s.variant)
    n = stash === nothing ? t.n : stash.nlive
    @inbounds for i in 1:n
        dmr = Int(t.dmr[i]); dmr == 0 && continue
        m = Float32(ie_dm_dg_mult(dgp, maxsp, Int(t.species[i]), dmr))
        m == 1f0 && continue
        t.diam_growth[i] *= m
        stash !== nothing && (stash.dgU[i] *= m; stash.dgL[i] *= m)
    end
    return
end

"""
    ie_dm_mortality_combine!(killed, s, fint, n)

mismrt.f:185-191: MAX-combine per-tree DM mortality (WKI = PROB·rate) into `killed[]`
(WK2 = max(WK2, WKI)) — DM mortality REPLACES background when larger, not additive.
No-op for non-IE / uninfected. Order-independent (per-tree max).
"""
function ie_dm_mortality_combine!(killed::AbstractVector{Float32}, s::StandState, fint::Float32, n::Int)
    _dm_effects_variant(s.variant) || return
    t = s.trees
    _, _, pmc, maxsp = _mis_tables(s.variant)
    @inbounds for i in 1:n
        dmr = Int(t.dmr[i]); dmr == 0 && continue
        pr = t.tpa[i]; pr <= 0f0 && continue
        rate = ie_dm_mortality_rate(pmc, maxsp, Int(t.species[i]), dmr, t.dbh[i], fint)
        wki = pr * rate
        killed[i] < wki && (killed[i] = wki)
    end
    return
end

"""
    ie_mistoe!(s; fint)

IE dwarf-mistletoe SPREAD/intensification for the cycle (mistoe.f MISTOE). Mutates `s.trees.dmr`
via the Hawksworth model: per host-species tree, prob-of-increase (PPLUS) intensifies existing
infection or introduces new infection into uninfected trees, prob-of-decrease (PMINUS) reduces it,
gated by the tallest-infected-tree height on the point (overstory vs understory branch). The
spread coefficients (BCONST/BDMR/BHTG/BTPA, DCONST/DDMR/DHTG/DTPA + the intensification thresholds)
are VARIANT-UNIFORM shared mistoe.f DATA (reused from the CR constants); only the host-species list
(IE_MIS_FIT) and species count (23) are IE-specific. YPLMLT/YNGMLT/DMMMLT default 1.0 (no MISTMULT).
Draws rann! per host tree ONLY when that species carries infection (SMR>0), matching FVS draw
count/order. No-op for non-IE. Mirrors the validated cr_mistoe! exactly.
"""
function ie_mistoe!(s::StandState; fint::Float32)
    _ie_mis_variant(s.variant) || return s
    t = s.trees
    t.n == 0 && return s
    species_sort!(s)
    isct = s.control.sp_count_tab
    ind1 = s.scratch.idx1
    rng  = s.rng
    fscale = fint / 10f0
    mis_fit, _, _, mis_maxsp = _mis_tables(s.variant)
    nsp = min(mis_maxsp, nspecies(s.variant))
    @inbounds for ispc in 1:nsp
        mis_fit[ispc] == 0 && continue
        i1 = isct[ispc, 1]; i1 == 0 && continue
        i2 = isct[ispc, 2]
        tottpa = 0f0; smr = 0f0
        for i3 in i1:i2
            i = Int(ind1[i3]); p = t.tpa[i]
            tottpa += p; smr += Float32(t.dmr[i]) * p
        end
        tottpa <= 0f0 && continue
        smr /= tottpa
        smr == 0f0 && continue                        # mistletoe-free species ⇒ NO draws
        dmtall = Dict{Int32,Float32}()
        for i3 in i1:i2
            i = Int(ind1[i3])
            if t.dmr[i] > 0
                pl = t.plot_id[i]; h = t.height[i]
                (get(dmtall, pl, 0f0) < h) && (dmtall[pl] = h)
            end
        end
        for i3 in i1:i2
            i = Int(ind1[i3])
            idmr = Int(t.dmr[i])
            htgr10 = t.ht_growth[i]
            pplus = 0f0
            if idmr < 6
                pplus = CR_DM_BCONST + CR_DM_BDMR[idmr + 1] +
                        CR_DM_BHTG * (htgr10 * 10f0 / fint) + CR_DM_BTPA * tottpa
                pplus != 0f0 && (pplus = 1f0 / (1f0 + fexp(-pplus)))
                pplus = pplus >= 1f0 ? 1f0 : 1f0 - fpow(1f0 - pplus, fscale)
            end
            if idmr != 0
                pminus = CR_DM_DCONST + CR_DM_DDMR * idmr +
                         CR_DM_DHTG * (htgr10 * 10f0 / fint) + CR_DM_DTPA * tottpa
                pminus != 0f0 && (pminus = 1f0 / (1f0 + fexp(-pminus)))
                pminus = pminus >= 1f0 ? 1f0 : 1f0 - fpow(1f0 - pminus, fscale)
                xnum = rann!(rng)
                dtall = get(dmtall, t.plot_id[i], 0f0)
                if idmr != 6 && pplus > xnum
                    if dtall * 0.7f0 > t.height[i]
                        x2 = rann!(rng)
                        m = Int(t.dmr[i])
                        inc = if m == 1
                            x2 < 0.61f0 ? 1 : (x2 < 0.83f0 ? 2 : 3)
                        elseif m == 2 || m == 3
                            x2 < 0.34f0 ? 1 : (x2 < 0.67f0 ? 2 : 3)
                        elseif m == 4
                            x2 < 0.55f0 ? 1 : 2
                        else
                            1
                        end
                        t.dmr[i] += Int32(inc)
                    else
                        t.dmr[i] += Int32(1)
                    end
                    t.dmr[i] > 6 && (t.dmr[i] = Int32(6))
                end
                (pminus > xnum) && (t.dmr[i] -= Int32(1))
            else
                xnum = rann!(rng)
                dtall = get(dmtall, t.plot_id[i], 0f0)
                if dtall * 0.7f0 > t.height[i]
                    if xnum < 0.55f0
                        x2 = rann!(rng)
                        t.dmr[i] = x2 < 0.69f0 ? Int32(1) : (x2 < 0.87f0 ? Int32(2) : Int32(3))
                    end
                else
                    (pplus > xnum) && (t.dmr[i] = Int32(1))
                end
            end
        end
        for i3 in i1:i2
            i = Int(ind1[i3])
            t.dmr[i] >= 4 && (t.mort_code[i] = Int32(3))
        end
    end
    return s
end
