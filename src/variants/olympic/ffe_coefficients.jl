# ffe_coefficients.jl (olympic) — OP FFE per-species coefficients (39 NWO ORGANON species).
# Extracted from bin/FVSop_buildDir: fmvinit.f (V2T), fmcrow.f (ISPMAP), fmcblk.f (BIOGRP=Jenkins),
# fmbrkt.f (bark), blkdat.f (ISPSPE sprout). Mirror of oregoncoast/ffe_coefficients.jl. sp38 = '---' placeholder.

const OP_FFE_V2T = Float32[24.9,23.1,21.8,19.3,22.5,20.6,23.1,26.2,21.8,20.6,23.7,21.2,21.2,22.5,23.7,28.1,21.2,19.3,26.2,26.2,27.4,23.1,36.2,36.2,36.2,21.8,19.3,37.4,34.9,29.9,22.5,23.7,26.2,27.4,27.4,29.3,22.5,28.1,29.9]
const OP_FFE_ISPMAP = Int32[4,4,4,1,4,18,4,8,20,18,11,15,15,15,13,3,19,7,6,24,5,23,10,17,17,41,17,17,16,1,14,11,7,56,57,61,64,0,41]
const OP_FFE_DBHMIN = 7.0f0
const OP_FFE_SPROUT = let v = zeros(Float32, 39); for sp in (17,21,22,23,24,25,26,27,28,33,34,35,36,37); v[sp] = 1f0; end; v end
