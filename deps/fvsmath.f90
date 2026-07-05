! fvsmath.f90 — FVSjl elementary-op Fortran companion (TOLERANCE_GOAL.md doctrine #8).
!
! The not-correctly-rounded Float32 transcendentals (exp, log, x**p) evaluate ~1 ULP
! differently between Julia's openlibm and gfortran's libm. FVS is built with gfortran,
! so to make FVSjl bit-identical to the live oracle in the compounding growth/crown/
! volume paths, Julia calls THESE routines (compiled by the same gfortran that builds
! bin/FVS*_buildDir) via ccall instead of its own openlibm. sqrt is IEEE-correctly-
! rounded in both and is intentionally NOT wrapped.
!
! Build (mirrors src/core/fmath.jl's _ensure_built):
!   gfortran -shared -fPIC -O2 -o libfvsmath.<dlext> fvsmath.f90
!
! Each function takes/returns REAL(4) BY VALUE with a C binding, matching the C ABI
! Julia's ccall expects.

function f32_exp(x) bind(C, name="f32_exp") result(y)
  use iso_c_binding, only: c_float
  real(c_float), value :: x
  real(c_float) :: y
  y = exp(x)
end function f32_exp

function f32_log(x) bind(C, name="f32_log") result(y)
  use iso_c_binding, only: c_float
  real(c_float), value :: x
  real(c_float) :: y
  y = log(x)
end function f32_log

function f32_pow(x, p) bind(C, name="f32_pow") result(y)
  use iso_c_binding, only: c_float
  real(c_float), value :: x, p
  real(c_float) :: y
  y = x ** p
end function f32_pow
