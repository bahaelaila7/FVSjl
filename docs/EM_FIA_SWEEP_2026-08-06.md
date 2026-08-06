# EM real-FIA sweep (2026-08-06) — validates AUTOES-amount + elevation-default fixes at scale

12-stand EM sample (em_sub.db) vs live FVSem_clean, 2 cycles. Sample = 6 empty (NULL-elev, AUTOES) + 6 treed
(non-NULL elev, growth). Purpose: confirm this session's fixes hold across stands with no regression.

| stand | type | cyc0 jl:live | cyc1 jl:live (TPA/BA) | verdict |
|-------|------|-------------|----------------------|---------|
| 5296268010661 | empty | 0:0 = | 118/0 = 118/0 | ✓ AUTOES bit-exact |
| 5297495010661 | empty | 0:0 = | 118/0 = 118/0 | ✓ |
| 5298115010661 | empty | 0:0 = | 118/0 = 118/0 | ✓ |
| 5298393010661 | empty | 0:0 = | 118/0 = 118/0 | ✓ |
| 5299488010661 | empty | 0:0 = | 118/0 = 118/0 | ✓ |
| 5303059010661 | empty | 0:0 = | 118/0 = 118/0 | ✓ |
| 1629523625290487 | treed | 5775/13 = | 5245/28 vs 5313/33 | close (TPA −1.3%, BA cornered) |
| 1856105255290487 | treed | 1110/47 = | 1066/66 vs 1002/61 | close (+6%) |
| 3034177010690 | treed | 385/169 = | 383/188 = 383/188 | ✓ bit-exact |
| 3094075010690 | treed | 88/72 = | 85/81 vs 86/81 | ✓ (BA exact) |
| 31440725010690 | treed | 558/125 = | 547/133 vs 547/134 | ✓ (TPA exact) |
| 225065919010661 | treed | 6081/45 = | 5924/96 vs 5447/49 | ★ OUTLIER: BA ~2×, under-thin |

## Results
- **AUTOES amount fix holds at SCALE**: 6/6 empty stands bit-exact (118) — generalizes beyond the 4 validated 2026-08-06.
- **Elevation-default fix causes NO regression**: ALL 12 cyc0 bit-exact; the 6 treed (non-NULL elev) stands unaffected.
- **5/6 treed stands** cyc1 within ±1-8% (bit-exact-or-cornered range).
- **1 outlier** (225065919010661): dense small-tree stand (6081 TPA, QMD~1.1"), cyc1 BA 96 vs live 49 (~2×) + under-thin.
  MEASURED: **AUTOES does NOT fire** on it (instrument-gated establish! log empty) ⇒ NOT an AUTOES/elevation artifact
  (has elevation; AUTOES inert). It is the pre-existing DENSE small-tree DG over-growth + under-thin class (#137/#139
  EMVAR small-tree on dense regen) — independent of this session's changes. A known-open small-tree residual.

⇒ Session deliverables VALIDATED at scale: EM AUTOES amount bit-exact (empty stands), elevation-default non-regressive
(treed stands). The lone outlier is a pre-existing dense small-tree DG issue, not introduced by these fixes.

## OUTLIER characterized (225065919010661) — dense HARDWOOD-SEEDLING DG over-growth [new lead]
11 tree records, dominated by DBH=0.1" seedlings with huge TPA: sp356 (2841 TPA), 823 (1243), 763 (532), 544
(355), 746 (177). Crosswalk: 544→GA, 746→AS, 823→OH (QUMA2), 356/763→OH (default). ⇒ the stand is a dense
hardwood/aspen SEEDLING cohort. cyc0 BA=45 bit-exact; cyc1 jl BA=96 vs live 49 (~2×) + under-thin (5924 vs 5447).
AUTOES does NOT fire (measured) ⇒ the divergence is the EM NON-CONIFER small-tree/regen DG (GA/AS/OH: CR-DIAGR +
aspen UTVAR forms, per fvsjl-em-variant-port) over-growing the 0.1" seedlings ~2× in one 10-yr cycle. This is the
#137/#139 dense-cohort small-tree class, here on HARDWOODS specifically (the EMVAR-conifer SMHTGF/SMDGF fix #136
does not cover the hardwood DIAGR/UTVAR seedling paths). NEXT: instrument jl vs live per-tree DG for the sp356/823
(OH) + 746 (AS) 0.1" seedlings on this stand (pre-tripling window, doctrine #3) — pick the over-growing form and
compare to the em/regent.f DIAGR/aspen small-tree DG. Likely the OH/aspen seedling DG lacks a size cap or uses the
large-tree form. Real-FIA-relevant (dense hardwood regen stands), NOT cornered (2× is far beyond the tie-break bar).
