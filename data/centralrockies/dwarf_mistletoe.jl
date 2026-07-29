# =============================================================================
# dwarf_mistletoe.jl (centralrockies) — CR dwarf mistletoe (DM) model data + pure kernels
#
# Ported from FVS mistoe/ (CR variant). DM is a WESTERN insect & pathogen (I&P)
# extension the eastern variants (SN/NE/CS/LS) never had, so the shared engine
# never needed it. FVS calls MISTOE unconditionally in the western base/gradd.f;
# FIA seeds per-tree DMR (Hawksworth 0-6 rating) from tree damage/severity codes,
# so DM fires automatically on infected stands with NO keyword.
#
# Discovered via the CR FIA sweep: the TPA-divergence class (jl over-retains) is
# substantially DM mortality — e.g. stand 46290437020004 (ponderosa, DMR 5.5,
# 91% infected) loses ~30%/decade to DM, which jl (lacking the model) never kills.
#
# Coefficients (misintcr.f) are IMODTY-INDEPENDENT for CR: MISFIT/DGPDMR/PMCSP are
# byte-identical across all 5 model-type blocks (verified). Growth+mortality coefs
# are noted "from the Utah variant" in misintcr.f. DMMMLT defaults 1.0, USEMRT=T.
#
# STATUS: data + pure kernels (RNG-INDEPENDENT) — validated where noted. The
# per-cycle spread/intensification (mistoe.f, RANN-heavy) + engine wiring
# (mortality apply, DGPDMR into diameter_growth!) are the remaining chunk; the
# spread's RNG stream-order is the bit-exact crux (shares seed 55329 / ZZRAN).
# =============================================================================

# --- MISFIT: which CR species (1..38) are DM hosts (misintcr.f AFIT, all types) ---
const CR_DM_MISFIT = Int8[
    1,1,1,1,1,1,0,1,1,1,   # 1-10   (7 RC not host)
    1,1,1,1,1,0,1,1,0,0,   # 11-20  (16 UJ, 19 WS, 20 AS not host)
    0,0,0,0,0,0,0,0,0,0,   # 21-30  hardwoods/junipers not host
    0,0,1,1,1,1,0,0]       # 31-38  (33 PM,34 PD,35 AZ,36 CI host)

# --- DGPDMR: diameter-growth potential multiplier by species × DMR(0..6) ---
# (misintcr.f ADGP; applied by misdgf: DG *= DGPDMR[sp][DMR+1]). Uniform over model types.
const CR_DM_DGPDMR = [
    (1.0,1.0,1.0,0.98,0.95,0.70,0.50), # 1  AF
    (1.0,1.0,1.0,0.98,0.95,0.70,0.50), # 2  CB
    (1.0,0.98,0.97,0.85,0.80,0.52,0.44), # 3  DF
    (1.0,1.0,1.0,0.98,0.95,0.70,0.50), # 4  GF
    (1.0,1.0,1.0,0.98,0.95,0.70,0.50), # 5  WF
    (1.0,1.0,1.0,0.98,0.86,0.73,0.50), # 6  MH
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 7  RC (not host)
    (1.0,0.94,0.92,0.88,0.84,0.58,0.54), # 8  WL
    (1.0,1.0,1.0,0.98,0.86,0.73,0.50), # 9  BC
    (1.0,1.0,1.0,1.0,0.94,0.80,0.59), # 10 LM
    (1.0,1.0,1.0,1.0,0.94,0.80,0.59), # 11 LP
    (1.0,1.0,1.0,1.0,0.94,0.80,0.59), # 12 PI
    (1.0,1.0,1.0,0.98,0.86,0.73,0.50), # 13 PP
    (1.0,1.0,1.0,1.0,0.94,0.80,0.59), # 14 WB
    (1.0,1.0,1.0,1.0,0.94,0.80,0.59), # 15 SW
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 16 UJ (not host)
    (1.0,1.0,1.0,0.98,0.86,0.73,0.50), # 17 BS
    (1.0,0.98,0.97,0.85,0.80,0.52,0.44), # 18 ES
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 19 WS (not host)
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 20 AS (not host)
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 21 NC
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 22 PW
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 23 GO
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 24 AW
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 25 EM
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 26 BK
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 27 SO
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 28 PB
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 29 AJ
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 30 RM
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 31 OJ
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 32 ER
    (1.0,1.0,1.0,1.0,0.94,0.80,0.59), # 33 PM
    (1.0,1.0,1.0,1.0,0.94,0.80,0.59), # 34 PD
    (1.0,1.0,1.0,1.0,0.94,0.80,0.59), # 35 AZ
    (1.0,1.0,1.0,0.98,0.86,0.73,0.50), # 36 CI
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0),     # 37 OS
    (1.0,1.0,1.0,1.0,1.0,1.0,1.0)]     # 38 OH

# --- PMCSP: DM mortality coefficients (B0,B1,B2) per species (misintcr.f APMC) ---
const CR_DM_PMCSP = [
    (0.0,0.00159,0.00508), (0.0,0.00159,0.00508), (0.01319,-0.01627,0.00822),
    (0.0,0.00159,0.00508), (0.0,0.00159,0.00508), (0.00681,-0.00580,0.00935),
    (0.0,0.0,0.0), (0.01319,-0.01627,0.00822), (0.00681,-0.00580,0.00935),
    (0.00112,0.02170,-0.00171), (0.00112,0.02170,-0.00171), (0.00112,0.02170,-0.00171),
    (0.00681,-0.00580,0.00935), (0.00112,0.02170,-0.00171), (0.00112,0.02170,-0.00171),
    (0.0,0.0,0.0), (0.00681,-0.00580,0.00935), (0.01319,-0.01627,0.00822),
    (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0),
    (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0),
    (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0), (0.0,0.0,0.0),
    (0.00112,0.02170,-0.00171), (0.00112,0.02170,-0.00171), (0.00112,0.02170,-0.00171),
    (0.00681,-0.00580,0.00935), (0.0,0.0,0.0), (0.0,0.0,0.0)]

# --- DMR spread probability coefficients (mistoe.f; same for all variants/species) ---
# FVS computes PPLUS/PMINUS in single-precision REAL, so these are Float32 (the spread is RNG-order sensitive).
const CR_DM_BDMR   = (0f0,2.45047f0,2.30723f0,1.8809f0,2.11457f0,1.43293f0,0f0)  # BDMR(DMR+1)
const CR_DM_BCONST = -1.67226f0;  const CR_DM_BTPA = -0.0012397f0; const CR_DM_BHTG = -0.0747205f0
const CR_DM_DDMR   =  0.0983757f0; const CR_DM_DCONST = -5.59798f0; const CR_DM_DTPA = -1.15053f-4
const CR_DM_DHTG   =  0.013267f0

# FIA Arceuthobium (dwarf mistletoe) damage codes (misdam.f). Legacy FVS codes 30-34
# carry the severity directly; these FIA codes map to severity (1-6) or default 3.
const CR_DM_FIA_CODES = Set{Int}([23005,23006,23007,23008,23009,23010,23011,23012,
    23013,23014,23015,23016,23017,23021,23023,23024])

"""
    cr_dm_seed_dmr(d1,s1,d2,s2,d3,s3) -> Int  (0..6)

Seed a tree's initial dwarf-mistletoe rating (IMIST) from its FIA damage/severity
code pairs (misdam.f). Legacy codes 30-34 → severity (clamped 0-6); FIA Arceuthobium
codes → severity if 1-6 else 3. First matching pair wins. VALIDATED bit-exact on
46290437020004 (10 records → DMR 6 from code 23023 sev 6, matching live's
"NUMBER OF RECORDS WITH MISTLETOE 10").
"""
function cr_dm_seed_dmr(d1,s1,d2,s2,d3,s3)
    @inbounds for (d,s) in ((d1,s1),(d2,s2),(d3,s3))
        if 30 <= d <= 34
            return clamp(Int(s),0,6)
        elseif d in CR_DM_FIA_CODES
            return (1 <= s <= 6) ? Int(s) : 3
        end
    end
    return 0
end

"""
    cr_dm_mortality_rate(sp, dmr, dbh, fint; dmmmlt=1.0) -> Float32

Periodic DM-induced mortality proportion for a tree (mismrt.f). Quadratic in DMR,
species multiplier, +20% for DBH<9, clamped, then annualized to the cycle length.
Added to the stand's per-tree mortality in the gradd apply. Returns 0 for DMR 0.
"""
function cr_dm_mortality_rate(sp::Integer, dmr::Integer, dbh::Real, fint::Real; dmmmlt::Real=1.0)
    dmr == 0 && return 0.0f0
    b0,b1,b2 = CR_DM_PMCSP[sp]
    m = b0 + b1*dmr + b2*dmr*dmr
    m *= dmmmlt
    small = dbh < 9.0
    small && (m *= 1.2)
    m < 0.0 && (m = 0.0)
    # DBH-dependent ceiling (mismrt.f:174-178): <9" caps at 0.71, >=9" caps at 0.5
    cap = small ? 0.71 : 0.5
    m > cap && (m = cap)
    return Float32(1.0 - (1.0 - m)^(fint/10.0))
end

"""
    cr_dm_dg_mult(sp, dmr) -> Float64   diameter-growth-potential multiplier (misdgf).
"""
cr_dm_dg_mult(sp::Integer, dmr::Integer) = CR_DM_DGPDMR[sp][dmr+1]
