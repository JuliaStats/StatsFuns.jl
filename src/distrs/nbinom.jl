# functions related to negative binomial distribution
# Implementation (log) pdf and cdfs and their inverses for the negative binomial distribution.
#
# Notation:
# - r > 0 is the number of successes until the experiment is stopped (generalized to real)
# - p ∈ [0,1] is the probability of a success (real)
# - k ≥ 0 is the number of failures (integer)
#
# The probability mass function is (k + r - 1 \choose k) p^r (1-p)^k

nbinompdf(r::Real, p::Real, k::Real) = exp(nbinomlogpdf(r, p, k))

nbinomlogpdf(r::Real, p::Real, k::Real) = nbinomlogpdf(promote(r, p, k)...)
function nbinomlogpdf(r::T, p::T, k::T) where {T <: Real}
    if !(0 <= p <= 1) || r <= 0
        return float(T)(NaN)
    elseif k >= 0 && isinteger(k)
        z = xlog1py(k, -p) + xlogy(r, p)
        iszero(k) && return z
        res = z - (logbeta(k + 1, r) + log(r + k))
        if !isnan(res)
            return res
        else
            return float(T)(-Inf)
        end
    else
        return float(T)(isnan(k) ? NaN : -Inf)
    end
end
function nbinomcdf(r::Real, p::Real, k::Real)
    if k < 0
        return zero(float(first(promote(r, p, k))))
    elseif isinf(k)
        return one(float(first(promote(r, p, k))))
    else
        return beta_inc(r, floor(k + 1), p)[1]
    end
end
function nbinomccdf(r::Real, p::Real, k::Real)
    if k < 0
        return one(float(first(promote(r, p, k))))
    elseif isinf(k)
        return zero(float(first(promote(r, p, k))))
    else
        return beta_inc(r, floor(k + 1), p)[2]
    end
end
function nbinomlogcdf(r::Real, p::Real, k::Real)
    if k < 0
        return oftype(float(first(promote(r, p, k))), -Inf)
    elseif isinf(k)
        return zero(float(first(promote(r, p, k))))
    else
        b1, b2 = beta_inc(r, floor(k + 1), p)
        return 10 * b1 < 7 ? log1p(-b2) : log(b1)
    end
end
function nbinomlogccdf(r::Real, p::Real, k::Real)
    if k < 0
        return zero(float(first(promote(r, p, k))))
    elseif isinf(k)
        return oftype(float(first(promote(r, p, k))), -Inf)
    else
        b1, b2 = beta_inc(r, floor(k + 1), p)
        return 10 * b1 < 7 ? log1p(-b1) : log(b2)
    end
end

# Inverse CDF: find smallest k such that nbinomcdf(r, p, k) >= cprob.
# Based on VBA critnegbinom by Ian Smith.
# PMF ratio: PMF(k+1)/PMF(k) = (k+r)/(k+1) * (1-p)
function _critnbinom(r::Float64, p::Float64, cprob::Float64)
    q = 1.0 - p

    # Normal approximation: mean = r*q/p, var = r*q/p^2
    μ = r * q / p
    σ = sqrt(μ / p)
    i = max(0.0, floor(μ + norminvcdf(min(cprob, 1.0 - 1.0e-15)) * σ + 0.5))

    # Compute CDF at the guess
    pr = Float64(nbinomcdf(r, p, i))

    if pr >= cprob
        # Search left
        while i > 0
            tpr = Float64(nbinompdf(r, p, i))
            if pr - tpr < cprob
                return i
            end
            pr -= tpr
            i -= 1.0
        end
        return 0.0
    else
        # Search right
        for _ in 1:10_000
            i += 1.0
            tpr = Float64(nbinompdf(r, p, i))
            pr += tpr
            if pr >= cprob
                return i
            end
        end
        return i
    end
end

# Inverse CCDF: find smallest k such that nbinomccdf(r, p, k) <= cprob.
# Based on VBA critcompnegbinom by Ian Smith.
function _critcompnbinom(r::Float64, p::Float64, cprob::Float64)
    q = 1.0 - p

    # Normal approximation
    μ = r * q / p
    σ = sqrt(μ / p)
    i = max(0.0, floor(μ - norminvcdf(min(cprob, 1.0 - 1.0e-15)) * σ + 0.5))

    # Compute CCDF at the guess
    pr = Float64(nbinomccdf(r, p, i))

    if pr > cprob
        # Search right
        for _ in 1:10_000
            i += 1.0
            tpr = Float64(nbinompdf(r, p, i))
            pr -= tpr
            if pr <= cprob
                return i
            end
        end
        return i
    else
        # Search left
        while i > 0
            tpr = Float64(nbinompdf(r, p, i))
            if pr + tpr > cprob
                return i
            end
            pr += tpr
            i -= 1.0
        end
        return 0.0
    end
end

# Wrapper with edge cases
function _nbinom_invcdf(r::Float64, p::Float64, q::Float64)
    if q < 0 || q > 1 || r <= 0 || p < 0 || p > 1 || isnan(q) || isnan(r) || isnan(p)
        return NaN
    elseif q == 0 || p == 1
        return 0.0
    elseif q == 1
        return Inf
    elseif p == 0
        return Inf
    end

    i = _critnbinom(r, p, q)

    # Post-correction
    pr = Float64(nbinomcdf(r, p, i))
    if pr >= q
        while i > 0
            pr2 = Float64(nbinomcdf(r, p, i - 1.0))
            if pr2 < q
                return i
            end
            i -= 1.0
        end
        return 0.0
    else
        return i + 1.0
    end
end

function _nbinom_invccdf(r::Float64, p::Float64, q::Float64)
    if q < 0 || q > 1 || r <= 0 || p < 0 || p > 1 || isnan(q) || isnan(r) || isnan(p)
        return NaN
    elseif q == 0 || p == 0
        return Inf
    elseif q == 1
        return 0.0
    elseif p == 1
        return 0.0
    end

    i = _critcompnbinom(r, p, q)

    # Post-correction
    pr = Float64(nbinomccdf(r, p, i))
    if pr <= q
        while i > 0
            pr2 = Float64(nbinomccdf(r, p, i - 1.0))
            if pr2 > q
                return i
            end
            i -= 1.0
        end
        return 0.0
    else
        return i + 1.0
    end
end

# Public API

function nbinominvcdf(r::Real, p::Real, q::Real)
    T = float(Base.promote_typeof(r, p, q))
    return convert(T, _nbinom_invcdf(Float64(r), Float64(p), Float64(q)))
end

function nbinominvccdf(r::Real, p::Real, q::Real)
    T = float(Base.promote_typeof(r, p, q))
    return convert(T, _nbinom_invccdf(Float64(r), Float64(p), Float64(q)))
end

function nbinominvlogcdf(r::Real, p::Real, lq::Real)
    T = float(Base.promote_typeof(r, p, lq))
    _lq = Float64(lq)
    result = if _lq > -1
        _nbinom_invccdf(Float64(r), Float64(p), -expm1(_lq))
    else
        _nbinom_invcdf(Float64(r), Float64(p), exp(_lq))
    end
    return convert(T, result)
end

function nbinominvlogccdf(r::Real, p::Real, lq::Real)
    T = float(Base.promote_typeof(r, p, lq))
    _lq = Float64(lq)
    result = if _lq > -1
        _nbinom_invcdf(Float64(r), Float64(p), -expm1(_lq))
    else
        _nbinom_invccdf(Float64(r), Float64(p), exp(_lq))
    end
    return convert(T, result)
end
