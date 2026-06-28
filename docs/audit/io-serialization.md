# Bandaid Audit — IO: .sum / .tre / DBS (+ stop-restart serialization)

Module files audited (FVSjl):
- `/workspace/FVSjl/src/io/summary.jl`  (SUMOUT `.sum` writer + row builder)
- `/workspace/FVSjl/src/io/dbs_output.jl` (DBS SQLite tables)
- `/workspace/FVSjl/src/io/treedata.jl` (INTREE FORMAT parse + `.tre` read/write)
- `/workspace/FVSjl/src/io/keyword.jl`  (KEYRDR lexer + `.key` write)
- `/workspace/FVSjl/src/io/input.jl`    (format-agnostic dispatch — pure routing, no FVS semantics)

FVS sources checked: `vbase/sumout.f`, `vbase/disply.f`, `vbase/evtstv.f`,
`base/keyrdr.f`, `sn/blkdat.f`, `fire/base/fmsdit.f`, `fire/base/fmmain.f`, the `dbs*.f` writers.

Scope note: the assignment lists "stop-restart serialization (putstd.f/getstd.f)", but there is
**no PUTSTD/GETSTD/stash serialization code in FVSjl at all** (`grep` for putstd/getstd/serialize/
fvsStopPoint over `src/` returns nothing in these files; the only `crown_lift`/snapshot machinery is
FFE state, not stop-restart). Stop/restart serialization is a FVSjulia (old port) feature; it is not
ported to FVSjl, so there is nothing here to flag for it. The `--stoppoint=` mechanism in FVS is a
command-line argument (base/cmdline.f:125), not a `.key` keyword.

---

## FLAG 1 — BANDAID: `snapshot_ffe_oldcrown!` seeded at inventory

- **jl symbol/line:** `write_sum_file`, `summary.jl:101-102`
  ```julia
  snapshot_ffe_oldcrown!(s)   # FMOLDC at inventory: gives the 1st cycle's crown-lift a valid
                              # OLDCRW (else the 1st cycle's fine down-wood is lost; DDW gap)
  ```
- **Claim:** without seeding the previous-cycle crown state (`ffe_old*`) at inventory, the *first*
  cycle's crown-lift fine down-wood would be lost (a DDW gap).
- **FVS source checked:** `fire/base/fmsdit.f:93` — `IF (ICYC.GT.1) THEN` guards the entire
  crown-lift / OLDCRW-falldown-to-down-wood block. fmsdit.f:64 also has `IF (ICYC .EQ. 1)` resetting
  fire state. FMSDIT is called from `base/grincr.f:227` once per cycle. So in **cycle 1 the crown-lift
  down-wood contribution is exactly zero by construction** — FVS deliberately produces no first-period
  crown-lift, because the "old crown" is undefined for the inventory records (the FVS comment at
  fmsdit.f:78-79 says so verbatim: "If we are on the first cycle, then the old crown is not known").
- **Why it's a bandaid:** seeding `ffe_old*` at inventory makes the first non-zero
  `compute_crown_lift!` (run at the end of the `c=0` iteration, `summary.jl:187`) compute a crown-base
  rise from the inventory→cycle-1 growth and book it as down wood in the following fuel loop — fabricating
  a crown-lift increment that FVS's `ICYC.GT.1` gate suppresses. The comment's justification ("else the
  1st cycle's down-wood is lost") is precisely the inverse of what fmsdit.f mandates: that down wood is
  *supposed* to be zero. This is the exact pattern flagged in the audit brief.
- **Faithfulness impact:** over-states DDW / fine-fuel loading by one period's worth of crown-lift in the
  first cycle of any FFE stand; perturbs FMCFMD small/large fuel-model selection and downstream fire
  behavior for stands that burn early. (The `fill!(crown_lift_annual,0)` at `:100` correctly zeroes the
  first *fuel loop*; the leak is the inventory `ffe_old*` seed feeding the first `compute_crown_lift!`.)

---

## FLAG 2 — GAP: `.sum` header sample-weight uses C `%E`, not Fortran `E15.7`

- **jl symbol/line:** `write_sum_header`, `summary.jl:58`
  ```julia
  @printf(io, "-999%5d %-26s %-4s%15.7E %-2s ...", ..., Float32(sample_wt), ...)
  ```
- **FVS source checked:** `vbase/sumout.f:107-109`
  ```fortran
  WRITE (JSUM2,2) LENG,NPLT,MGMID,SAMWT,VARACD,DAT,TIM,...
2 FORMAT ('-999',I5,1X,A26,1X,A4,E15.7,5(1X,A),I3)
  ```
  FVS writes SAMWT with the Fortran `E15.7` edit descriptor, which normalizes the mantissa to `0.x`
  (e.g. `1.0 → "  0.1000000E+01"`, 7 significant digits). Julia's `@printf "%15.7E"` follows C and emits
  `d.dddddddE±dd` (e.g. `1.0 → "  1.0000000E+00"`; verified live: `0.04 → "  3.9999999E-02"`). The two
  forms differ in **both** the exponent and the digit before the point for every value.
- **Severity:** GAP — the period data rows are bit-exact (`_SUM_ROW_FMT` matches sumout.f format `9014`,
  verified field-by-field), and the `-999` header is a marker line that the `.sum` row comparisons in the
  harness do not byte-check, so this is latent. But a faithful `.sum` header would not match FVS for any
  stand whose `SamplingWt ≠` a value that happens to coincide — i.e. effectively always.
- **Faithfulness impact:** the SAMWT token of the per-stand `-999` header diverges from FVS byte-for-byte;
  any consumer that parses the header sample weight as fixed-format Fortran-E text would misread it.

---

## FLAG 3 — GAP (latent): `STOP` detected on a 4-char prefix, not the full 8-col field

- **jl symbol/line:** `read_keyword!`, `keyword.jl:117`
  ```julia
  head8[1:4] == "STOP" && return _record("STOP", record, ..., KW_STOP, 0)
  ```
- **FVS source checked:** `base/keyrdr.f:60` — `TMP=RECORD(1:8)` (CHARACTER*8), then
  `IF (TMP.EQ.'STOP')`. Fortran pads the literal `'STOP'` to the operand width (8), so the test is
  `RECORD(1:8) .EQ. 'STOP    '` — it fires **only** when columns 1-8 are exactly `STOP` followed by four
  blanks. FVSjl fires on any record whose first four columns are `STOP`, e.g. `STOPXXXX` or a hypothetical
  `STOPPNT`/`STOPPOINT` keyword.
- **Severity:** GAP, currently **latent** — a scan of `base/*.f` and `sn/*.f` keyword literals finds no
  FVS keyword that begins with `STOP` other than `STOP` itself (the stop-point feature is the
  `--stoppoint=` command-line arg, not a card), so no real keyword collides today. It is still a genuine
  deviation from keyrdr.f:60's exact-field semantics and would mis-terminate a keyword stream if any
  `STOP`-prefixed card (or a free-form supplemental line starting with `STOP`) ever appeared.
- **Faithfulness impact:** none on the current SN keyword set; a correctness trap if the keyword
  vocabulary grows. Tighten to require cols 1-8 == `"STOP    "` (i.e. `head8 == "STOP    "`).

---

## FLAG 4 — GAP: `FVS_Carbon` removed / fire-released columns hard-zeroed

- **jl symbol/line:** `write_dbs_carbon!`, `dbs_output.jl:129`
  ```julia
  ... Float64(r.total), 0.0, 0.0))   # Total_Removed_Carbon, Carbon_Released_From_Fire
  ```
- **Claim (comment, :114-115):** "Total-Removed / Released-from-Fire are 0 (no harvest / fire carbon
  accounting on the carbon-report path yet)."
- **FVS source checked:** `dbsfmcrpt.f` (the DBSFMCRPT writer feeds these two columns from the FFE carbon
  accounting, non-zero whenever the cycle has harvest removals or a fire event). FVSjl always writes 0.
- **Severity:** GAP — faithful for the no-harvest / no-fire carbon scenario actually exercised, but a
  carbon report on a stand WITH a thin or a SIMFIRE burn would under-report Total_Removed_Carbon and
  Carbon_Released_From_Fire as 0 versus FVS's computed values. The comment honestly admits the
  incompleteness; flagged so it is tracked as a real divergence rather than a parity field.
- **Faithfulness impact:** two `FVS_Carbon` columns wrong (always 0) for managed / burned stands.

---

## Items reviewed and found FAITHFUL (not individually flagged)

- `_SUM_ROW_FMT` (`summary.jl:16-19`) — matches sumout.f:222 FORMAT `9014`
  (`2I4,I6,I4,I5,2I4,F5.1,9I6,I4,I5,2I4,F5.1,2X,I6,I5,I6,2X,F6.1,1X,I3,1X,2I1`) field-for-field.
- MAI terminal-row quirk (`summary.jl:117-121, 230-233`) — the cited lines verify: evtstv.f:414
  `IF(AGE.GT.AGELST)TOTREM=IOSUM(9,ICYC-1)+TOTREM` (one-cycle lag; SN falls to CASE DEFAULT, the
  `'SN'` branch at :387-393 is commented out), disply.f:392
  `BCYMAI=(IOSUM(5,IKNT)+TOTREM)/IOSUM(2,IKNT)`. The "final row excludes the last growing cycle's
  removal" is the genuine consequence of the lag; the jl `cum_rem_merch - prev_increment` on the last
  row reproduces it.
- `summary_row` NINT rounding `trunc(Int, x+0.5)` — matches FVS `INT(.../GROSPC+0.5)` (disply.f:355-359).
- RESETAGE rebasing (`summary.jl:226-228`) — consistent with resage.f age-rebase semantics.
- `DEFAULT_TREE_FORMAT` (`treedata.jl:19-20`) — byte-identical to `sn/blkdat.f:58-59` `TREFMT`.
- The `I4,T1,I7` plot/id overlap and implied-decimal F parsing — intentional FVS behavior, reproduced.
- KEYRDR control flow (`!` skip, pre-heading blank skip, `*`/blank skip after heading, COMMENT…END
  consume) — matches keyrdr.f:50-100. `_scan_parms` NF = `div(pos-11,10)` matches keyrdr.f
  `IP=K+IP-11; NF=(IP-1)/10`. CR/CRLF/CR-only normalization is a reasonable cross-platform read of what
  Fortran's per-OS record reader handles.
- DBS schema column orders / NINT (`trunc(x+0.5f0)`) for SDImax & site index, and the dynamic
  `FVS_Compute` table — consistent with the cited `dbs*.f` writers. `FVS_Cases` build metadata is
  environment-specific (not a parity field) as documented.
- `write_tree_file` / `write_keyfile` overlap-anchor heuristics — these are FVSjl's *own* round-trip
  writers (FVS has no `.tre`/`.key` writer); a semantic round-trip, not a claim of FVS-output parity, so
  not a faithfulness concern.
