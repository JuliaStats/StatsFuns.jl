# functions related to student's T distribution

tdistpdf(ν::Real, x::Real) = exp(tdistlogpdf(ν, x))

tdistlogpdf(ν::Real, x::Real) = tdistlogpdf(promote(ν, x)...)
function tdistlogpdf(ν::T, x::T) where {T <: Real}
    isinf(ν) && return normlogpdf(x)
    νp12 = (ν + 1) / 2
    return loggamma(νp12) - (logπ + log(ν)) / 2 - loggamma(ν / 2) - νp12 * log1p(x^2 / ν)
end

function tdistcdf(ν::T, x::T) where {T <: Real}
    if isinf(ν)
        return normcdf(x)
    elseif x < 0
        return fdistccdf(one(ν), ν, x^2) / 2
    else
        return 1 - fdistccdf(one(ν), ν, x^2) / 2
    end
end
tdistcdf(ν::Real, x::Real) = tdistcdf(map(float, promote(ν, x))...)

tdistccdf(ν::Real, x::Real) = tdistcdf(ν, -x)

function tdistlogcdf(ν::T, x::T) where {T <: Real}
    if isinf(ν)
        return normlogcdf(x)
    elseif x < 0
        ret = fdistlogccdf(one(ν), ν, x^2)
        return ret - logtwo
    else
        return log1p(-fdistccdf(one(ν), ν, x^2) / 2)
    end
end
tdistlogcdf(ν::Real, x::Real) = tdistlogcdf(map(float, promote(ν, x))...)

tdistlogccdf(ν::Real, x::Real) = tdistlogcdf(ν, -x)

# Pure Julia implementation of the inverse CDF, based on VBA invtdist by Ian Smith:
# start from a Cornish-Fisher expansion of the quantile around the normal quantile
# (central region) or from inverting the leading power-law term of the tail (tail
# region), then polish with Newton's method applied to log(cdf), where each
# iteration costs a single evaluation of the incomplete beta function.

# Crossover curves between the Cornish-Fisher and tail starting approximations.
# Based on VBA BetterThanTailApprox by Ian Smith.
function _tdist_cf_start_is_better(pr::Float64, ν::Float64)
    if ν <= 2.0
        return pr > 0.25 * exp((1.0 - ν) * 1.78514841051368)
    elseif ν <= 5.0
        return pr > 0.045 * exp((2.0 - ν) * 1.30400766847605)
    elseif ν <= 20.0
        return pr > 0.0009 * exp((5.0 - ν) * 0.921034037197618)
    else
        return pr > 9.0e-10 * exp((20.0 - ν) * 0.690775527898214)
    end
end

# cdf and pdf of the t distribution at x <= 0, with B = beta(ν / 2, 1 / 2) precomputed
function _tdistcdf_pdf(ν::Float64, B::Float64, x::Float64)
    if abs(x) >= min(1.0, ν)
        # this form of k2 = ν / (ν + x^2) and x2 = x^2 / (ν + x^2) avoids
        # premature overflow of x^2
        k2 = ν / x
        t = x + k2
        k2 = k2 / t
        x2 = x / t
    else
        x² = x * x
        t = ν + x²
        x2 = x² / t
        k2 = ν / t
    end
    p = first(beta_inc(ν / 2, 0.5, k2, x2)) / 2
    f = k2^((ν + 1) / 2) / (sqrt(ν) * B)
    return p, f
end

function _tdistinvcdf(ν::Float64, p::Float64)
    if isnan(ν) || isnan(p) || !(0.0 <= p <= 1.0) || !(ν > 0.0)
        return NaN
    elseif isinf(ν)
        return norminvcdf(p)
    elseif p == 0.0
        return -Inf
    elseif p == 1.0
        return Inf
    elseif p == 0.5
        return 0.0
    end

    # work with the smaller tail; the result is negated for p > 1/2 on return
    pr = p > 0.5 ? 1.0 - p : p
    small = 1.0e-14
    smalllpr = -small * log(pr) * pr

    logB = logbeta(ν / 2, 0.5)
    B = exp(logB)

    local tp::Float64
    tprob = 0.0
    tpdif = 0.0
    if pr >= 0.5 || (ν >= 1.0 && _tdist_cf_start_is_better(pr, ν))
        # Cornish-Fisher expansion of the t quantile in powers of 1/ν around the
        # normal quantile (Fisher & Cornish, Technometrics 2 (1960), 209-225)
        xn = norminvcdf(pr)
        x = xn * xn
        tp = (((((27.0 * x + 339.0) * x + 930.0) * x - 1782.0) * x - 765.0) * x + 17955.0) / (368640.0 * ν)
        tp = (tp + ((((79.0 * x + 776.0) * x + 1482.0) * x - 1920.0) * x - 945.0) / 92160.0) / ν
        tp = (tp + (((3.0 * x + 19.0) * x + 17.0) * x - 15.0) / 384.0) / ν
        tp = (tp + ((5.0 * x + 16.0) * x + 3.0) / 96.0) / ν
        tp = (tp + (x + 1.0) / 4.0) / ν
        tp = xn * (1.0 + tp)
        # for large ν the expansion has already converged to full precision and the
        # Newton polish can be skipped (validated relative error <= 2.2e-15 in this
        # region against a high-precision reference)
        if ν >= 250.0 && x <= ν / 100.0
            return p > 0.5 ? -tp : tp
        end
        tprob = 0.0
        tpdif = 1.0 + abs(tp)
    elseif ν < 1.0
        # leading power-law term of the tail solved in log space:
        # pr ≈ (ν / t^2)^(ν / 2) / (ν * B)
        lν = log(ν)
        tp = -exp(lν / 2 - (log(pr) + lν + logB) / ν)
        isfinite(tp) || return p > 0.5 ? -tp : tp
        tprob, f = _tdistcdf_pdf(ν, B, tp)
        if f < floatmin(Float64)
            tpdif = 0.0
        else
            tpdif = tprob / f * log1p((tprob - pr) / pr)
            tpnew = tp - tpdif
            tp = tpnew < 0.0 ? tpnew : tp / 2
        end
    else
        # invert the leading power-law term of the tail: pr ≈ f(0) * √ν * |t|^(-ν)
        # where f(0) = 1 / (√ν * B) is the density at zero, map back to the t scale,
        # apply one closed-form next-order correction, and one Newton step
        u = exp(-log(ν * B * pr) / ν)
        tp = -sqrt(ν) * sqrt(u - 1.0) * sqrt(u + 1.0)
        isfinite(tp) || return p > 0.5 ? -tp : tp
        tpdif = tp / ν
        tpdif = -log1p((0.5 - 1.0 / (ν + 2.0)) / (1.0 + tpdif * tp)) * (tpdif + 1.0 / tp)
        tp -= tpdif
        tprob, f = _tdistcdf_pdf(ν, B, tp)
        if f < floatmin(Float64)
            tpdif = 0.0
        else
            tpdif = tprob / f * log1p((tprob - pr) / pr)
            tpnew = tp - tpdif
            tp = tpnew < 0.0 ? tpnew : tp / 2
        end
    end

    # Newton iteration applied to log(cdf): step = cdf / pdf * log(cdf / pr).
    # Near the center this is an ordinary Newton step; in the tails, where the cdf
    # is nearly exponential in t, the log transform makes the problem nearly linear.
    iter = 0
    while abs(tprob - pr) > smalllpr && abs(tpdif) > small * (1.0 + abs(tp))
        (iter += 1) > 100 && break
        tprob, f = _tdistcdf_pdf(ν, B, tp)
        f < floatmin(Float64) && break
        tpdif = tprob / f * log1p((tprob - pr) / pr)
        tpnew = tp - tpdif
        # keep the iterate in the negative half-line, where the solution lies
        tp = tpnew < 0.0 ? tpnew : tp / 2
    end
    return p > 0.5 ? -tp : tp
end

# The kernel operates in Float64, like the Rmath-based functions elsewhere in the
# package; argument types with more precision than Float64 are truncated.
function tdistinvcdf(ν::T, p::T) where {T <: Real}
    return convert(float(T), _tdistinvcdf(Float64(ν), Float64(p)))
end
tdistinvcdf(ν::Real, p::Real) = tdistinvcdf(promote(ν, p)...)

tdistinvccdf(ν::Real, p::Real) = -tdistinvcdf(ν, p)
function tdistinvlogcdf(ν::T, logp::T) where {T <: Real}
    if isinf(ν)
        return norminvlogcdf(logp)
    else
        logq = logp + logtwo
        if logq < 0
            return -sqrt(fdistinvlogccdf(one(ν), ν, logq))
        else
            return sqrt(fdistinvlogccdf(one(ν), ν, log2mexp(logq)))
        end
    end
end
tdistinvlogcdf(ν::Real, logp::Real) = tdistinvlogcdf(map(float, promote(ν, logp))...)

tdistinvlogccdf(ν::Real, logp::Real) = -tdistinvlogcdf(ν, logp)
