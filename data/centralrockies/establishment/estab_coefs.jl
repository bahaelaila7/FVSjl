# CR establishment coefficients (cr/blkdat.f). Per-species (38, in the CR species order).
# XMIN = establishment min height (estab_min_ht); HHTMAX = max established height cap (esgent.f).
# Ready for the establish!(::CentralRockies) branch: the CR regen/planted trees take XMIN initial
# height and grow via the CR REGENT model (esgent.f uses REGENT — already ported as _cr_regent_tree),
# capped at HHTMAX — NOT the eastern ESSUBH height-at-age path in establish!.
const _CR_ES_XMIN = Float32[
    0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5,
    1.0, 0.5, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 3.0,
    3.0, 3.0, 0.5, 0.5, 0.5, 0.5, 0.5, 3.0, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.5]
const _CR_ES_HHTMAX = Float32[
    7.0, 7.0, 10.0, 7.0, 7.0, 10.0, 9.0, 10.0, 9.0, 9.0,
    10.0, 6.0, 10.0, 9.0, 9.0, 6.0, 7.0, 7.0, 7.0, 16.0,
    16.0, 16.0, 10.0, 10.0, 10.0, 10.0, 10.0, 16.0, 6.0, 6.0,
    6.0, 6.0, 6.0, 6.0, 6.0, 10.0, 9.0, 12.0]
