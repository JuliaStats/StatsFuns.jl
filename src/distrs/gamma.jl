# functions related to gamma distribution

using HypergeometricFunctions: drummond1F1

# R implementations
using .RFunctions:
    # gammapdf,
    # gammalogpdf,
    # gammacdf,
    # gammaccdf,
    # gammalogcdf,
    # gammalogccdf,
    # gammainvcdf,
    # gammainvccdf,
    gammainvlogcdf,
    gammainvlogccdf

# Julia implementations
gammapdf(k::Real, θ::Real, x::Real) = exp(gammalogpdf(k, θ, x))

gammalogpdf(k::Real, θ::Real, x::Real) = gammalogpdf(promote(k, θ, x)...)
function gammalogpdf(k::T, θ::T, x::T) where {T<:Real}
    # we ensure that `log(x)` does not error if `x < 0`
    xθ = max(x, 0) / θ
    val = -loggamma(k) + xlogy(k - 1, xθ) - log(θ) - xθ
    return x < 0 ? oftype(val, -Inf) : val
end

gammacdf(k::T, θ::T, x::T) where {T<:Real} = first(gamma_inc(k, x/θ))
gammacdf(k::Real, θ::Real, x::Real)        = gammacdf(promote(float(k), θ, x)...)

gammaccdf(k::T, θ::T, x::T) where {T<:Real} = last(gamma_inc(k, x/θ))
gammaccdf(k::Real, θ::Real, x::Real)        = gammaccdf(promote(float(k), θ, x)...)

gammalogcdf(k::Real, θ::Real, x::Real) = _gammalogcdf(map(float, promote(k, θ, x))...)

# Implemented via the non-log version. For tiny values, we recompute the result with
# loggamma. In that situation, there is little risk of significant cancellation.
function _gammalogcdf(k::Float64, θ::Float64, x::Float64)
    l, u = gamma_inc(k, x/θ)
    if l < eps(Float64)
        return -log(k) + k*log(x/θ) - x/θ + log(drummond1F1(1.0, 1 + k, x/θ)) - loggamma(k)
    elseif l < 0.7
        return log(l)
    else
        return log1p(-u)
    end
end
function _gammalogcdf(k::T, θ::T, x::T) where {T<:Union{Float16,Float32}}
    return T(_gammalogcdf(Float64(k), Float64(θ), Float64(x)))
end

gammalogccdf(k::Real, θ::Real, x::Real) = _gammalogccdf(map(float, promote(k, θ, x))...)

# Implemented via the non-log version. For tiny values, we recompute the result with
# loggamma. In that situation, there is little risk of significant cancellation.
function _gammalogccdf(k::Float64, θ::Float64, x::Float64)
    l, u = gamma_inc(k, x/θ)
    if u < eps(Float64)
        return loggamma(k, x/θ) - loggamma(k)
    elseif u < 0.7
        return log(u)
    else
        return log1p(-l)
    end
end

function _gammalogccdf(k::T, θ::T, x::T) where {T<:Union{Float16,Float32}}
    return T(_gammalogccdf(Float64(k), Float64(θ), Float64(x)))
end

function gammainvcdf(k::Real, θ::Real, p::Real)
    _k, _θ, _p = map(float, promote(k, θ, p))
    return _θ*gamma_inc_inv(_k, _p, 1 - _p)
end

function gammainvccdf(k::Real, θ::Real, p::Real)
    _k, _θ, _p = map(float, promote(k, θ, p))
    return _θ*gamma_inc_inv(_k, 1 - _p, _p)
end
