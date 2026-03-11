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

# Rmath implementations
# TODO: implement https://arxiv.org/abs/2001.03953
# for inverting the incomplete beta function wrt the 2nd argument.
function nbinominvcdf(r::Real, p::Real, q::Real)
    T = float(Base.promote_typeof(r, p, q))
    return convert(T, Rmath.qnbinom(q, r, p, true, false))
end
function nbinominvccdf(r::Real, p::Real, q::Real)
    T = float(Base.promote_typeof(r, p, q))
    return convert(T, Rmath.qnbinom(q, r, p, false, false))
end
function nbinominvlogcdf(r::Real, p::Real, lq::Real)
    T = float(Base.promote_typeof(r, p, lq))
    return convert(T, Rmath.qnbinom(lq, r, p, true, true))
end
function nbinominvlogccdf(r::Real, p::Real, lq::Real)
    T = float(Base.promote_typeof(r, p, lq))
    return convert(T, Rmath.qnbinom(lq, r, p, false, true))
end
