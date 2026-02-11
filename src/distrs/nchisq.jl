# functions related to noncentral chi-square distribution

# R implementations
using .RFunctions:
    nchisqinvcdf,
    nchisqinvccdf,
    nchisqinvlogcdf,
    nchisqinvlogccdf

nchisqpdf(k::Real, λ::Real, x::Real) = nchisqpdf(promote(k, λ, x)...)
function nchisqpdf(k::T, λ::T, x::T) where T<:Real
    half = one(T) / 2
    rootxλ = sqrt(x * λ)
    bix = besselix(k / 2 - 1, rootxλ)
    return exp(rootxλ - x / 2 - λ / 2 + xlogy(k / 4 - half, x / λ)) * bix / 2
end

nchisqlogpdf(k::Real, λ::Real, x::Real) = nchisqlogpdf(promote(k, λ, x)...)
function nchisqlogpdf(k::T, λ::T, x::T) where T<:Real
    half = one(T) / 2
    rootxλ = sqrt(x * λ)
    logbix = log(besselix(k / 2 - 1, rootxλ))
    return log(half) + rootxλ - x / 2 - λ / 2 + xlogy(k / 4 - half, x / λ) + logbix
end

_nchisqcdf(k::Real, λ::Real, x::Real, invert::Bool) = _nchisqcdf(promote(k, λ, x)..., invert)
function _nchisqcdf(k::T, λ::T, x::T, invert::Bool) where T<:Real
    s = u = exp(-λ / 2)
    t = exp(log(x / 2) * k / 2 - x / 2 - loggamma(k / 2 + 1))
    res = s * t
    i = 1
    if x < k + λ # cdf < ~0.5 so don't invert
        while s * t > res * eps(T)
            u *= λ / (2i)
            s += u
            t *= x / (k + 2i)
            res = muladd(s, t, res)
            i += 1
        end
    else
        invert = !invert
        res = fma(-s, t, 1)
        while s * t > res * eps(T)
            u *= λ / (2i)   
            s += u
            t *= x / (k + 2i)
            res = fma(-s, t, res)
            i += 1
        end
    end
    return res, invert
end

function nchisqcdf(k::Real, λ::Real, x::Real)
    res, invert = _nchisqcdf(k, λ, x, false)
    invert ? 1 - res : res
end
function nchisqccdf(k::Real, λ::Real, x::Real)
    res, invert = _nchisqcdf(k, λ, x, true)
    invert ? 1 - res : res
end
function nchisqlogcdf(k::Real, λ::Real, x::Real)
    res, invert = _nchisqcdf(k, λ, x, false)
    invert ? log1p(-res) : log(res)
end

function nchisqlogccdf(k::Real, λ::Real, x::Real)
    res, invert = _nchisqcdf(k, λ, x, true)
    invert ? log1p(-res) : log(res)
end
