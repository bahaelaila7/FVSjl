# Divergences from Fortran FVS

This file records **every** place where FVSjl knowingly differs from the reference
Fortran FVS, in two categories:

1. **Suspected FVS bugs** — behaviour we believe is wrong given the code's intent or
   comments. These are fixed only behind the non-default `faithful=false` path and
   listed here so they can be raised with the FVS maintainers for confirmation.
2. **Unavoidable numerical drift** — e.g. a transcendental function whose last ulp
   differs between the Julia and Fortran math libraries. Documented, not chased.

In `faithful=true` mode (the default, and what the test suite runs) FVSjl is
intended to be bit-exact to Fortran, so this list has no behavioural effect there.

Format per entry:

> ### <short title>
> - **Where (FVSjl):** file:func
> - **Where (Fortran):** file.f:line
> - **Category:** suspected-bug | numeric-drift
> - **Gate:** `faithful=false` (bug fixes only) | always (drift)
> - **Description / evidence:**
> - **Status:** open | reported | confirmed

---

_No divergences recorded yet. Entries are added as porting uncovers them._
