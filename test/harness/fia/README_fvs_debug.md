# Rebuilding a debug-enabled FVSls binary (breaks the per-tree-volume tooling wall)

The relinked /tmp/FVSls_new can't recompile individual .f (gfortran .mod version mismatch).
Workaround — recompile the small module chain + the target file, link a SEPARATE binary:

    BD=/workspace/ForestVegetationSimulator/bin/FVSls_buildDir
    mkdir -p /tmp/dbgbuild && cd /tmp/dbgbuild
    cp "$BD"/debug_mod.f "$BD"/clkcoef_mod.f "$BD"/r9clark.f "$BD"/r9coeff.inc .
    # add WRITE(0,*) debug lines to r9clark.f (e.g. in r9bdft / after `call r9bdft` where dbhOb is in scope)
    gfortran -c -O2 -fno-automatic debug_mod.f clkcoef_mod.f
    gfortran -c -O2 -fno-automatic r9clark.f
    OBJS=$(ls "$BD"/*.o | grep -vE "/(debug_mod|clkcoef_mod|r9clark)\.o$")
    gfortran -o /tmp/FVSls_dbg $OBJS debug_mod.o clkcoef_mod.o r9clark.o /tmp/glibc_shim.o -lpthread -ldl

/tmp/FVSls_dbg reproduces the oracle .sum exactly; debug goes to stderr (unit 0). This is how the
LS-hardwood BdFt residual was localized (r9bdft per-log board, r9clark per-tree BDTREE dbh+board).

## Extending the debug binary to the volume interface (fvsvol.f / volinit.f)

The r9clark recipe generalizes. To dump the board-type dispatch (Slice S43 — the
`VOL(2)=VOL(10)` International swap), recompile `fvsvol.f` (+ `mrules_mod.f`) — it only
needs the `.F77` COMMON includes from the buildDir, no module-version cascade:

```
BD=/workspace/ForestVegetationSimulator/bin/FVSls_buildDir
cd /tmp/fia_val/dbgbuild
cp "$BD"/fvsvol.f "$BD"/mrules_mod.f "$BD"/*.F77 .
# add e.g.  WRITE(0,*)'NATDBG',ISPC,METHB(ISPC),D,BBFV,TVOL(2),TVOL(10)  after the
# 'IF(BBFV.LT.0.)BBFV=0.' line in the NATCRS entry (keep the line ≤ col 72, fixed-form!)
gfortran -c -O2 -fno-automatic mrules_mod.f && gfortran -c -O2 -fno-automatic fvsvol.f
OBJS=$(ls "$BD"/*.o | grep -vE "/(fvsvol|mrules_mod)\.o$")
gfortran -o /tmp/FVSls_dbg2 $OBJS fvsvol.o mrules_mod.o /tmp/glibc_shim.o -lpthread -ldl
```

`volinit.f` (where the actual `VOL(2)=VOL(10)` swap lives) pulls in `charmod`→`volinput_mod`→…
and IS a version cascade — recompile the whole chain, or (easier) confirm the swap indirectly
from `fvsvol.f`'s NATCRS dump: if `TVOL(2)==TVOL(10)` on arrival, volinit already swapped it.
Fixed-form gotcha: keep every added `WRITE` within column 72 or it silently truncates and the
build serves the STALE binary (bit us once — the `,VOLEQ` tail split into a phantom `VO` symbol).
