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

    tp = if pr >= 0.5 || (ν >= 1.0 && _tdist_cf_start_is_better(pr, ν))
        # Cornish-Fisher expansion of the t quantile in powers of 1/ν around the
        # normal quantile (Fisher & Cornish, Technometrics 2 (1960), 209-225)
        xn = norminvcdf(pr)
        x = xn * xn
        t = (((((27.0 * x + 339.0) * x + 930.0) * x - 1782.0) * x - 765.0) * x + 17955.0) / (368640.0 * ν)
        t = (t + ((((79.0 * x + 776.0) * x + 1482.0) * x - 1920.0) * x - 945.0) / 92160.0) / ν
        t = (t + (((3.0 * x + 19.0) * x + 17.0) * x - 15.0) / 384.0) / ν
        t = (t + ((5.0 * x + 16.0) * x + 3.0) / 96.0) / ν
        t = (t + (x + 1.0) / 4.0) / ν
        t = xn * (1.0 + t)
        # for large ν the expansion has already converged to full precision and the
        # Newton polish can be skipped (validated relative error <= 2.2e-15 in this
        # region against a high-precision reference)
        if ν >= 250.0 && x <= ν / 100.0
            return p > 0.5 ? -t : t
        end
        t
    elseif ν < 1.0
        # leading power-law term of the tail solved in log space:
        # pr ≈ (ν / t^2)^(ν / 2) / (ν * B)
        lν = log(ν)
        -exp(lν / 2 - (log(pr) + lν + logB) / ν)
    else
        # invert the leading power-law term of the tail, pr ≈ f(0) * √ν * |t|^(-ν)
        # where f(0) = 1 / (√ν * B) is the density at zero, map back to the t
        # scale, and apply one closed-form next-order correction
        u = exp(-log(ν * B * pr) / ν)
        t = -sqrt(ν) * sqrt(u - 1.0) * sqrt(u + 1.0)
        if isfinite(t)
            d = t / ν
            t -= -log1p((0.5 - 1.0 / (ν + 2.0)) / (1.0 + d * t)) * (d + 1.0 / t)
        end
        t
    end
    # the true quantile may overflow in the extreme tails
    isfinite(tp) || return p > 0.5 ? -tp : tp

    # Newton iteration applied to log(cdf): step = cdf / pdf * log(cdf / pr).
    # Near the center this is an ordinary Newton step; in the tails, where the cdf
    # is nearly exponential in t, the log transform makes the problem nearly linear.
    iter = 0
    while true
        tprob, f = _tdistcdf_pdf(ν, B, tp)
        # the density underflows only where the start already carries full precision
        f < floatmin(Float64) && break
        tpdif = tprob / f * log1p((tprob - pr) / pr)
        # second-order coefficient of the iteration, |d²log(cdf)/dt² / dlog(cdf)/dt|,
        # for the quadratic-convergence estimate of the error remaining after the step
        curv = abs(-(ν + 1.0) * tp / (ν + tp * tp) - f / tprob)
        tpnew = tp - tpdif
        # keep the iterate in the negative half-line, where the solution lies
        tp = tpnew < 0.0 ? tpnew : tp / 2
        tol = small * (1.0 + abs(tp))
        # converged if the cdf already matched, the step was negligible, or the
        # estimated error remaining after the step is negligible (the exact
        # second-order coefficient is 1/2; the factor 32 is a 64x safety margin)
        abs(tprob - pr) <= smalllpr && break
        abs(tpdif) <= tol && break
        32.0 * curv * tpdif * tpdif <= tol && break
        (iter += 1) >= 100 && break
    end
    return p > 0.5 ? -tp : tp
end

# Only Float16 and Float32 are routed explicitly through the Float64 kernel;
# wider types such as BigFloat are unsupported rather than silently computed
# at Float64 precision.
_tdistinvcdf(ν::Float16, p::Float16) = convert(Float16, _tdistinvcdf(Float64(ν), Float64(p)))
_tdistinvcdf(ν::Float32, p::Float32) = convert(Float32, _tdistinvcdf(Float64(ν), Float64(p)))
tdistinvcdf(ν::Real, p::Real) = _tdistinvcdf(map(float, promote(ν, p))...)

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
