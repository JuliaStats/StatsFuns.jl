# functions related to gamma distribution

using HypergeometricFunctions: _₁F₁

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
    val = -loggamma(k) - log(θ) - xθ
    # xlogy(k - 1, xθ) - xθ -> -∞ for xθ -> ∞ so we only add the first term
    # when it's safe
    if isfinite(xθ)
        val += xlogy(k - 1, xθ)
    end
    return x < 0 ? oftype(val, -Inf) : val
end

function gammacdf(k::T, θ::T, x::T) where {T<:Real}
    # Handle the degenerate case
    if iszero(k)
        return float(last(promote(x, x >= 0)))
    end
    return first(gamma_inc(k, max(0, x)/θ))
end
gammacdf(k::Real, θ::Real, x::Real) = gammacdf(map(float, promote(k, θ, x))...)

function gammaccdf(k::T, θ::T, x::T) where {T<:Real}
    # Handle the degenerate case
    if iszero(k)
        return float(last(promote(x, x < 0)))
    end
    return last(gamma_inc(k, max(0, x)/θ))
end
gammaccdf(k::Real, θ::Real, x::Real) = gammaccdf(map(float, promote(k, θ, x))...)

gammalogcdf(k::Real, θ::Real, x::Real) = _gammalogcdf(map(float, promote(k, θ, x))...)

# Implemented via the non-log version. For tiny values, we recompute the result with
# loggamma. In that situation, there is little risk of significant cancellation.
function _gammalogcdf(k::Float64, θ::Float64, x::Float64)
    # Handle the degenerate case
    if iszero(k)
        return log(x >= 0)
    end

    xdθ = max(0, x)/θ
    l, u = gamma_inc(k, xdθ)
    if l < floatmin(Float64) && isfinite(k) && isfinite(xdθ)
        return -log(k) + k*log(xdθ) - xdθ + log(_₁F₁(1.0, 1 + k, xdθ)) - loggamma(k)
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
    # Handle the degenerate case
    if iszero(k)
        return log(x < 0)
    end

    xdθ = max(0, x)/θ
    l, u = gamma_inc(k, xdθ)
    if u < floatmin(Float64)
        return loggamma(k, xdθ) - loggamma(k)
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
    _k, _θ, _p = promote(k, θ, p)
    return _θ*gamma_inc_inv(_k, _p, 1 - _p)
end

function gammainvccdf(k::Real, θ::Real, p::Real)
    _k, _θ, _p = promote(k, θ, p)
    return _θ*gamma_inc_inv(_k, 1 - _p, _p)
end
