# Functions related to hypergeometric distribution
# Pure Julia implementation based on VBA code by Ian Smith
# https://iandjmsmith.wordpress.com/
# License: MIT

# Constants
const _hyper_cfVSmall = 1.0e-15
const _hyper_scalefactor = 1.1579208923731619542357098500869e+77   # 2^256
const _hyper_scalefactor2 = 8.6361685550944446253863518628004e-78  # 2^-256
const _hyper_minLog1Value = -0.79149064
const _hyper_max_discrete = 9.007199254740991e15  # 2^53
const _hyper_max_crit = 4.503599627370496e15      # 2^52

# Internal PMF computation
# ai = x, aji = n - x, aki = ms - x, amkji = mf - n + x
function _hypergeometric_term(ai::Float64, aji::Float64, aki::Float64, amkji::Float64)
    ak = aki + ai       # ms
    amk = amkji + aji   # mf
    aj = aji + ai       # n
    am = amk + ak       # ms + mf
    amj = amkji + aki   # ms + mf - n

    if am > _hyper_max_discrete
        return NaN
    end

    if ai == 0 && (aji <= 0 || aki <= 0 || amj < 0 || amk < 0)
        return 1.0
    elseif ai > 0 && min(aki, aji) == 0 && max(amj, amk) == 0
        return 1.0
    elseif ai >= 0 && amkji > -1 && aki > -1 && aji >= 0
        c1 = lfbaccdif1(ak, amk) - lfbaccdif1(ai, aki) - lfbaccdif1(ai, aji) - lfbaccdif1(aki, amkji) - logfbit(ai)

        ai1 = ai + 1.0; aj1 = aj + 1.0; ak1 = ak + 1.0; am1 = am + 1.0
        aki1 = aki + 1.0; aji1 = aji + 1.0
        amk1 = amk + 1.0; amj1 = amj + 1.0; amkji1 = amkji + 1.0

        cjkmi = ab_minus_cd(aji, aki, ai, amkji)

        c5 = (cjkmi - ai) / (amkji1 * am1)
        c3 = if c5 < _hyper_minLog1Value
            amkji * (log((amj1 * amk1) / (amkji1 * am1)) - c5) - c5
        else
            amkji * log1pmx(c5) - c5
        end

        c5 = (-cjkmi - aji) / (aki1 * am1)
        c4 = if c5 < _hyper_minLog1Value
            aki * (log((ak1 * amj1) / (aki1 * am1)) - c5) - c5
        else
            aki * log1pmx(c5) - c5
        end
        c3 += c4

        c5 = (-cjkmi - aki) / (aji1 * am1)
        c4 = if c5 < _hyper_minLog1Value
            aji * (log((aj1 * amk1) / (aji1 * am1)) - c5) - c5
        else
            aji * log1pmx(c5) - c5
        end
        c3 += c4

        c5 = (cjkmi - amkji) / (ai1 * am1)
        c4 = if c5 < _hyper_minLog1Value
            ai * (log((aj1 * ak1) / (ai1 * am1)) - c5) - c5
        else
            ai * log1pmx(c5) - c5
        end
        c3 += c4

        logterm = (c1 + 1.0 / am1) + c3
        sqrtterm = sqrt((amk1 * ak1) * (aj1 * amj1) / ((amkji1 * aki1 * aji1) * (am1 * ai1)))
        return exp(logterm) * sqrtterm * Float64(invsqrt2π)
    else
        return 0.0
    end
end

# Internal CDF computation
function _hypergeometric(ai::Float64, aji::Float64, aki::Float64, amkji::Float64, comp::Bool)
    # Determine swap direction for numerical stability
    if amkji > -1 && amkji < 0
        ip1 = -amkji
        mkji = ip1 - 1.0
        allIntegral = false
    else
        ip1 = amkji + 1.0
        mkji = amkji
        allIntegral = ai == floor(ai) && aji == floor(aji) && aki == floor(aki) && mkji == floor(mkji)
    end

    if allIntegral
        swapped = (ai + 0.5) * (mkji + 0.5) >= (aki - 0.5) * (aji - 0.5)
    elseif (ai < 100 && ai == floor(ai)) || mkji < 0
        swapped = if comp
            (ai + 0.5) * (mkji + 0.5) >= aki * aji
        else
            (ai + 0.5) * (mkji + 0.5) >= aki * aji + 1000
        end
    elseif ai < 1
        swapped = (ai + 0.5) * (mkji + 0.5) >= aki * aji
    elseif aji < 1 || aki < 1 || (ai < 1 && ai > 0)
        swapped = false
    else
        swapped = (ai + 0.5) * (mkji + 0.5) >= (aki - 0.5) * (aji - 0.5)
    end

    if !swapped
        i = ai; ji = aji; ki = aki
    else
        i = aji - 1.0; ji = ai + 1.0; ki = ip1
        ip1 = aki; mkji = aki - 1.0
    end

    c2 = ji + i
    c4_pop = mkji + ki + c2  # population size

    if c4_pop > _hyper_max_discrete
        return NaN
    end

    if (i >= 0 && (ji <= 0 || ki <= 0)) || (ip1 + ki <= 0) || (ip1 + ji <= 0)
        exact = true
        prob = i >= 0 ? 1.0 : 0.0
    elseif ip1 > 0 && ip1 < 1
        exact = false
        prob = _hypergeometric_term(i, ji, ki, ip1) * (ip1 * (c4_pop + 1.0)) / ((ki + ip1) * (ji + ip1))
    else
        exact = (i == 0 && (ji <= 0 || ki <= 0 || mkji + ki < 0 || mkji + ji < 0)) ||
            (i > 0 && min(ki, ji) == 0 && max(mkji + ki, mkji + ji) == 0)
        prob = _hypergeometric_term(i, ji, ki, mkji)
    end

    if exact || prob == 0.0
        return (swapped == comp) ? prob : 1.0 - prob
    end

    a1 = 0.0
    c4 = c4_pop  # working copy for CF

    if i < mkji
        must_do_cf = i != floor(i)
        maxSums = floor(i)
    else
        must_do_cf = mkji != floor(mkji)
        maxSums = floor(max(mkji, 0.0))
    end

    if must_do_cf
        sumAlways = 0
        sumFactor = 5
    else
        sumAlways = 20
        sumFactor = 10
    end

    if maxSums > sumAlways || must_do_cf
        numb = floor(sumFactor / c4 * exp(log((ki + i) * (ji + i) * (ip1 + ji) * (ip1 + ki)) / 3.0))
        numb = floor(i - (ki + i) * (ji + i) / c4 + numb)
        numb = clamp(numb, 0.0, maxSums)
    else
        numb = maxSums
    end

    if 2 * numb <= maxSums || must_do_cf
        # Continued fraction evaluation
        b1 = 1.0
        c1 = 0.0
        c2_cf = i - numb
        c3 = mkji - numb
        s = c3
        a2 = c2_cf
        c3 -= 1.0
        b2 = ab_minus_cd(ki + numb + 1.0, ji + numb + 1.0, c2_cf - 1.0, c3)
        bn = b2
        bnAdd = c3 + c4 + c2_cf - 2.0

        while b2 > 0 && abs(a2 * b1 - a1 * b2) > abs(_hyper_cfVSmall * b1 * a2)
            c1 += 1.0; c2_cf -= 1.0
            an = (c1 * c2_cf) * (c3 * c4)
            c3 -= 1.0; c4 -= 1.0
            bn += bnAdd; bnAdd -= 4.0
            a1 = bn * a2 + an * a1
            b1 = bn * b2 + an * b1
            if b1 > _hyper_scalefactor
                a1 *= _hyper_scalefactor2; b1 *= _hyper_scalefactor2
                a2 *= _hyper_scalefactor2; b2 *= _hyper_scalefactor2
            end
            c1 += 1.0; c2_cf -= 1.0
            an = (c1 * c2_cf) * (c3 * c4)
            c3 -= 1.0; c4 -= 1.0
            bn += bnAdd; bnAdd -= 4.0
            a2 = bn * a1 + an * a2
            b2 = bn * b1 + an * b2
            if b2 > _hyper_scalefactor
                a1 *= _hyper_scalefactor2; b1 *= _hyper_scalefactor2
                a2 *= _hyper_scalefactor2; b2 *= _hyper_scalefactor2
            end
        end

        if b1 < 0 || b2 < 0
            return NaN
        else
            a1 = a2 / b2 * s
        end
    else
        numb = maxSums
    end

    # Direct summation
    c1_s = i - numb + 1.0
    c2_s = mkji - numb + 1.0
    c3_s = ki + numb
    c4_s = ji + numb
    for _ in 1:Int(numb)
        a1 = (1.0 + a1) * ((c1_s * c2_s) / (c3_s * c4_s))
        c1_s += 1.0; c2_s += 1.0; c3_s -= 1.0; c4_s -= 1.0
    end

    a1 = (1.0 + a1) * prob

    if swapped == comp
        return a1
    else
        return a1 > 0.99 ? NaN : 1.0 - a1
    end
end

# Inverse CDF search (lower tail)
# Faithful port of VBA crithyperg by Ian Smith
const _hyper_nearly_zero = 9.99999983659714e-317
const _hyper_cSmall = 5.562684646268003457725581793331e-309

function _crithyperg(j::Float64, k::Float64, m::Float64, cprob::Float64)
    if cprob > 0.5
        return _critcomphyperg(j, k, m, 1.0 - cprob)
    end

    mx = min(j, k)
    mn = max(0.0, j + k - m)

    # Normal approximation for initial guess
    i = j * k / m + norminvcdf(cprob) * sqrt(j * k * (m - j) * (m - k) / (m * m * max(m - 1.0, 1.0)))

    while true
        i = clamp(floor(i + 0.5), mn, mx)
        if i >= _hyper_max_crit
            return i
        end
        pr = _hypergeometric(i, j - i, k - i, m - k - j + i, false)
        tpr = 0.0
        if pr >= cprob
            if i == mn
                return mn
            end
            tpr = _hypergeometric_term(i, j - i, k - i, m - k - j + i)
            if pr < 1.00001 * tpr
                # PMF dominates: ratio stepping left
                tpr *= ((i + 1.0) * (m - j - k + i + 1.0)) / ((k - i) * (j - i))
                i -= 1.0
                while tpr > cprob
                    tpr *= ((i + 1.0) * (m - j - k + i + 1.0)) / ((k - i) * (j - i))
                    i -= 1.0
                end
                # Falls through to top of while(true) for re-evaluation
            else
                pr -= tpr
                if pr < cprob
                    return i
                end
                i -= 1.0
                if i == mn
                    return mn
                end
                temp = (pr - cprob) / tpr
                if temp > 10
                    # Large jump with parabolic refinement
                    temp = floor(temp + 0.5)
                    i -= temp
                    temp2 = _hypergeometric_term(i, j - i, k - i, m - k - j + i)
                    i -= temp * (tpr - temp2) / (2.0 * temp2)
                else
                    tpr *= ((i + 1.0) * (m - j - k + i + 1.0)) / ((k - i) * (j - i))
                    pr -= tpr
                    if pr < cprob
                        return i
                    end
                    i -= 1.0
                    temp2 = (pr - cprob) / tpr
                    if temp2 < temp - 0.9
                        # Linear stepping
                        while pr >= cprob
                            tpr *= ((i + 1.0) * (m - j - k + i + 1.0)) / ((k - i) * (j - i))
                            pr -= tpr
                            i -= 1.0
                        end
                        return i + 1.0
                    else
                        # Log-ratio jump
                        ratio = ((i + 1.0) * (m - j - k + i + 1.0)) / ((k - i) * (j - i))
                        temp = floor(log(cprob / pr) / log(ratio) + 0.5)
                        i -= temp
                        temp2 = _hypergeometric_term(i, j - i, k - i, m - k - j + i)
                        if temp2 > _hyper_nearly_zero
                            temp = log((cprob / pr) * (tpr / temp2)) / log(ratio)
                            i -= temp
                        end
                    end
                end
                # Falls through to top of while(true) for re-evaluation
            end
        else
            # Search right
            while tpr < _hyper_cSmall && pr < cprob
                i += 1.0
                tpr = _hypergeometric_term(i, j - i, k - i, m - k - j + i)
                pr += tpr
            end
            while pr < cprob
                i += 1.0
                tpr *= ((k - i + 1.0) * (j - i + 1.0)) / (i * (m - j - k + i))
                pr += tpr
            end
            return i
        end
    end
    return
end

# Inverse CDF search (upper tail)
# Faithful port of VBA critcomphyperg by Ian Smith
function _critcomphyperg(j::Float64, k::Float64, m::Float64, cprob::Float64)
    if cprob > 0.5
        return _crithyperg(j, k, m, 1.0 - cprob)
    end

    mx = min(j, k)
    mn = max(0.0, j + k - m)

    # Normal approximation
    i = j * k / m - norminvcdf(cprob) * sqrt(j * k * (m - j) * (m - k) / (m * m * max(m - 1.0, 1.0)))

    while true
        i = clamp(floor(i + 0.5), mn, mx)
        if i >= _hyper_max_crit
            return i
        end
        pr = _hypergeometric(i, j - i, k - i, m - k - j + i, true)
        tpr = 0.0
        if pr > cprob
            i += 1.0
            tpr = _hypergeometric_term(i, j - i, k - i, m - k - j + i)
            if pr < (1.0 + 0.00001) * tpr
                # PMF dominates: ratio stepping right
                while tpr > cprob
                    i += 1.0
                    tpr *= ((k - i + 1.0) * (j - i + 1.0)) / (i * (m - j - k + i))
                end
                # Falls through to top of while(true) for re-evaluation
            else
                pr -= tpr
                if pr <= cprob
                    return i
                end
                temp = (pr - cprob) / tpr
                if temp > 10
                    # Large jump with parabolic refinement
                    temp = floor(temp + 0.5)
                    i += temp
                    temp2 = _hypergeometric_term(i, j - i, k - i, m - k - j + i)
                    i += temp * (tpr - temp2) / (2.0 * temp2)
                else
                    i += 1.0
                    tpr *= ((k - i + 1.0) * (j - i + 1.0)) / (i * (m - j - k + i))
                    pr -= tpr
                    if pr <= cprob
                        return i
                    end
                    temp2 = (pr - cprob) / tpr
                    if temp2 < temp - 0.9
                        # Linear stepping
                        while pr > cprob
                            i += 1.0
                            tpr *= ((k - i + 1.0) * (j - i + 1.0)) / (i * (m - j - k + i))
                            pr -= tpr
                        end
                        return i
                    else
                        # Log-ratio jump
                        ratio = ((k - i + 1.0) * (j - i + 1.0)) / (i * (m - j - k + i))
                        temp = floor(log(cprob / pr) / log(ratio) + 0.5)
                        i += temp
                        temp2 = _hypergeometric_term(i, j - i, k - i, m - k - j + i)
                        if temp2 > _hyper_nearly_zero
                            temp = log((cprob / pr) * (tpr / temp2)) / log(ratio)
                            i += temp
                        end
                    end
                end
                # Falls through to top of while(true) for re-evaluation
            end
        else
            # Search left
            while tpr < _hyper_cSmall && pr <= cprob
                tpr = _hypergeometric_term(i, j - i, k - i, m - k - j + i)
                pr += tpr
                i -= 1.0
            end
            while pr <= cprob
                tpr *= ((i + 1.0) * (m - j - k + i + 1.0)) / ((k - i) * (j - i))
                pr += tpr
                i -= 1.0
            end
            return i + 1.0
        end
    end
    return
end

# Wrappers matching VBA crit_hypergeometric / comp_crit_hypergeometric
function _hyper_invcdf(ms::Float64, mf::Float64, n::Float64, p::Float64)
    m = ms + mf
    mn = max(0.0, n - mf)
    mx = min(ms, n)

    if p < 0 || p > 1 || isnan(p)
        return NaN
    elseif p == 0
        return mn
    elseif ms == 0 || n == 0
        return 0.0
    elseif mf == 0 || n == m
        return mx
    elseif p == 1
        return mx
    end

    i = _crithyperg(n, ms, m, p)

    # Post-correction (from crit_hypergeometric)
    pr = _hypergeometric(i, n - i, ms - i, mf - n + i, false)
    if pr == p
        return i
    elseif pr > p
        i2 = i - 1.0
        if i2 >= mn
            pr2 = _hypergeometric(i2, n - i2, ms - i2, mf - n + i2, false)
            if pr2 >= p
                return i2
            end
        end
        return i
    else
        return i + 1.0
    end
end

function _hyper_invccdf(ms::Float64, mf::Float64, n::Float64, q::Float64)
    m = ms + mf
    mn = max(0.0, n - mf)
    mx = min(ms, n)

    if q < 0 || q > 1 || isnan(q)
        return NaN
    elseif q == 1
        return mn
    elseif ms == 0 || n == 0
        return 0.0
    elseif mf == 0 || n == m
        return mx
    elseif q == 0
        return mx
    end

    i = _critcomphyperg(n, ms, m, q)

    # Post-correction (from comp_crit_hypergeometric)
    pr = _hypergeometric(i, n - i, ms - i, mf - n + i, true)
    if pr == q
        return i
    elseif pr < q
        i2 = i - 1.0
        if i2 >= mn
            pr2 = _hypergeometric(i2, n - i2, ms - i2, mf - n + i2, true)
            if pr2 <= q
                return i2
            end
        end
        return i
    else
        return i + 1.0
    end
end

# Public API

function hyperpdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    _x = round(Float64(x))
    _ms = Float64(ms); _mf = Float64(mf); _n = Float64(n)
    result = _hypergeometric_term(_x, _n - _x, _ms - _x, _mf - _n + _x)
    return convert(T, result)
end

function hyperlogpdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    return convert(T, log(Float64(hyperpdf(ms, mf, n, x))))
end

function hypercdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    _x = floor(Float64(x) + 1.0e-7)
    _ms = Float64(ms); _mf = Float64(mf); _n = Float64(n)
    result = _hypergeometric(_x, _n - _x, _ms - _x, _mf - _n + _x, false)
    return convert(T, result)
end

function hyperccdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    _x = floor(Float64(x) + 1.0e-7)
    _ms = Float64(ms); _mf = Float64(mf); _n = Float64(n)
    result = _hypergeometric(_x, _n - _x, _ms - _x, _mf - _n + _x, true)
    return convert(T, result)
end

function hyperlogcdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    return convert(T, log(Float64(hypercdf(ms, mf, n, x))))
end

function hyperlogccdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    return convert(T, log(Float64(hyperccdf(ms, mf, n, x))))
end

function hyperinvcdf(ms::Real, mf::Real, n::Real, q::Real)
    T = float(Base.promote_typeof(ms, mf, n, q))
    result = _hyper_invcdf(Float64(ms), Float64(mf), Float64(n), Float64(q))
    return convert(T, result)
end

function hyperinvccdf(ms::Real, mf::Real, n::Real, q::Real)
    T = float(Base.promote_typeof(ms, mf, n, q))
    result = _hyper_invccdf(Float64(ms), Float64(mf), Float64(n), Float64(q))
    return convert(T, result)
end

function hyperinvlogcdf(ms::Real, mf::Real, n::Real, lq::Real)
    T = float(Base.promote_typeof(ms, mf, n, lq))
    _lq = Float64(lq)
    isinf(_lq) && return convert(T, NaN)
    result = _hyper_invcdf(Float64(ms), Float64(mf), Float64(n), exp(_lq))
    return convert(T, result)
end

function hyperinvlogccdf(ms::Real, mf::Real, n::Real, lq::Real)
    T = float(Base.promote_typeof(ms, mf, n, lq))
    _lq = Float64(lq)
    isinf(_lq) && return convert(T, NaN)
    result = _hyper_invccdf(Float64(ms), Float64(mf), Float64(n), exp(_lq))
    return convert(T, result)
end
