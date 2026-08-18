const EPS  = 1.0e-6
const BMIN = 1.0e-10

h2d(s) = reinterpret(Float64, parse(UInt64, s, base=16))

# MPBGAM (mpbgam.f): double gamma via Stirling log-gamma (shift arg >=18).
function mpbgam(xx::Float64)::Float64
    zz = xx
    local dlgam::Float64
    if xx > 1.0e10
        dlgam = xx < 1.0e70 ? zz*(log(zz)-1.0) : 1.0e75
    elseif xx <= 1.0e-9
        dlgam = -1.0e75
    else
        term = 1.0
        while zz <= 18.0
            term *= zz; zz += 1.0
        end
        rz2 = 1.0/zz^2
        dlgam = (zz-0.5)*log(zz) - zz + Float64(0.9189385332046727f0) - log(term) +
            (1.0/zz)*(0.8333333333333333e-1 - (rz2*(0.2777777777777777e-2 +
            (rz2*(0.7936507936507936e-3 - (rz2*(0.5952380952380952e-3)))))))
    end
    return exp(dlgam)
end

# PQSML (pqsml.f): power-series for the incomplete beta.
function pqsml(x::Float64, p::Float64, q::Float64)::Float64
    u = x^p; s = u/p; xk = 0.0
    while true
        xk += 1.0
        u = (xk-q)*x*u/xk
        v = u/(xk+p)
        s += v
        abs(v) > abs(EPS*s) || break
    end
    return s*mpbgam(p+q)/(mpbgam(p)*mpbgam(q))
end

# FORW (forw.f): forward recurrence.
function forw(x::Float64, p::Float64, q::Float64, x1::Float64, x2::Float64, n::Int)::Float64
    n < 1 && return x1
    n == 1 && return x2
    n1 = n-1; t1 = x1; pq1 = p+q-1.0; xa1 = 1.0-x; forr = x2
    for i in 1:n1
        xi = Float64(i); temp = forr
        forr = forr + (forr-t1)*(pq1+xi)*xa1/(xi+q)
        t1 = temp
    end
    return forr
end

# BACK (back.f): backward continued fraction.
function back(x::Float64, p::Float64, q::Float64, z::Float64, n::Int)::Float64
    bk = z
    n == 0 && return bk
    q1 = q-1.0; ap = z; fn = Float64(n); xnu = Float64(n+n+6)
    while true
        xn = xnu; r = 0.0; bk = z
        while true
            xn -= 1.0
            t = p + xn
            u = (t+q1)*x
            r = u/(u+t-t*r)
            xn <= fn && (bk = r*bk)
            bk < BMIN && return bk
            xn > 1.0 || break
        end
        abs(bk-ap) < abs(EPS*bk) && return bk
        ap = bk; xnu += 5.0
    end
end

# BETIN (betin.f): main dispatcher.
function betin(a::Float64, b::Float64, x::Float64)::Float64
    (a > 0.0 && b > 0.0 && x >= 0.0 && x <= 1.0) || return 0.0
    ia = trunc(Int, a); p = a - ia
    if p == 0.0; p = 1.0; ia -= 1; end
    (x == 0.0 || x == 1.0) && return x
    local bet::Float64
    if x <= 0.5
        q = b; m = trunc(Int, q); q1 = q - m
        if q1 == 0.0; q1 = 1.0; m -= 1; end
        x1 = pqsml(x, p, q1)
        x2 = m > 0 ? pqsml(x, p, q1+1.0) : 0.0
        bet = back(x, p, q, forw(x, p, q1, x1, x2, m), ia)
    else
        y = 1.0 - x; q = p; p = b; m = trunc(Int, p); p1 = p - m
        if p1 == 0.0; m -= 1; p1 = 1.0; end
        x1 = back(y, p1, q,      pqsml(y, p1, q),      m)
        x2 = back(y, p1, q+1.0,  pqsml(y, p1, q+1.0),  m)
        bet = 1.0 - forw(y, p, q, x1, x2, ia)
    end
    bet < BMIN && (bet = 0.0)
    return bet
end

