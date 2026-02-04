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
    return loggamma(r + x) - loggamma(x + 1) - loggamma(r) + log1p(-p) * x + xlogy(r, p)
end
function nbinomcdf(r::Real, p::Real, x::Real)
    if x < 0
        return zero(float(first(promote(r, p, x))))
    elseif isinf(x)
        return one(float(first(promote(r, p, x))))
    end
    return beta_inc(r, floor(x + 1), p)[1]
end
function nbinomccdf(r::Real, p::Real, x::Real)
    if x < 0
        return one(float(first(promote(r, p, x))))
    elseif isinf(x)
        return zero(float(first(promote(r, p, x))))
    end
    return beta_inc(r, floor(x + 1), p)[2]
end
function nbinomlogcdf(r::Real, p::Real, x::Real)
    if x < 0
        return oftype(float(first(promote(r, p, x))), -Inf)
    elseif isinf(x)
        return zero(float(first(promote(r, p, x))))
    end
    b1, b2 = beta_inc(r, floor(x + 1), p)
    return 10 * b1 < 7 ? log1p(-b2) : log(b1)
end
function nbinomlogccdf(r::Real, p::Real, x::Real)
    if x < 0
        return zero(float(first(promote(r, p, x))))
    elseif isinf(x)
        return oftype(float(first(promote(r, p, x))), -Inf)
    end
    b1, b2 = beta_inc(r, floor(x + 1), p)
    return 10 * b1 < 7 ? log1p(-b1) : log(b2)
end

# TODO: implement https://arxiv.org/abs/2001.03953
# for inverting the incomplete beta function wrt the 2nd argument.
