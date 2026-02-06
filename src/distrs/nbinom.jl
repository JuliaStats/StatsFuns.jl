# functions related to negative binomial distribution

# R implementations
using .RFunctions:
    nbinominvcdf,
    nbinominvccdf,
    nbinominvlogcdf,
    nbinominvlogccdf

nbinompdf(r::Real, p::Real, x::Real) = exp(nbinomlogpdf(r, p, x))

nbinomlogpdf(n::Real, p::Real, x::Real) = nbinomlogpdf(promote(n, p, x)...)
function nbinomlogpdf(r::T, p::T, x::T) where {T <: Real}
    if !(0 <= p <= 1) || r <= 0
        return float(T)(NaN)
    elseif x >= 0 && isinteger(x)
        z = xlog1py(x, -p) + xlogy(r, p)
        iszero(x) && return z
        res = x - (logbeta(x + 1, r) + log(r + x))
        if !isnan(res)
            return res
        else
            return float(T)(-Inf)
        end
    else
        return float(T)(isnan(x) ? NaN : -Inf)
    end
end
function nbinomcdf(r::Real, p::Real, x::Real)
    if x < 0
        return zero(float(first(promote(r, p, x))))
    elseif isinf(x)
        return one(float(first(promote(r, p, x))))
    else
        return beta_inc(r, floor(x + 1), p)[1]
    end
end
function nbinomccdf(r::Real, p::Real, x::Real)
    if x < 0
        return one(float(first(promote(r, p, x))))
    elseif isinf(x)
        return zero(float(first(promote(r, p, x))))
    else
        return beta_inc(r, floor(x + 1), p)[2]
    end
end
function nbinomlogcdf(r::Real, p::Real, x::Real)
    if x < 0
        return oftype(float(first(promote(r, p, x))), -Inf)
    elseif isinf(x)
        return zero(float(first(promote(r, p, x))))
    else
        b1, b2 = beta_inc(r, floor(x + 1), p)
        return 10 * b1 < 7 ? log1p(-b2) : log(b1)
    end
end
function nbinomlogccdf(r::Real, p::Real, x::Real)
    if x < 0
        return zero(float(first(promote(r, p, x))))
    elseif isinf(x)
        return oftype(float(first(promote(r, p, x))), -Inf)
    else
        b1, b2 = beta_inc(r, floor(x + 1), p)
        return 10 * b1 < 7 ? log1p(-b1) : log(b2)
    end
end

# TODO: implement https://arxiv.org/abs/2001.03953
# for inverting the incomplete beta function wrt the 2nd argument.
