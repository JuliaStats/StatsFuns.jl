# functions related to gamma distribution

using HypergeometricFunctions: _₁F₁

# Julia implementations
gammapdf(k::Real, θ::Real, x::Real) = exp(gammalogpdf(k, θ, x))

gammalogpdf(k::Real, θ::Real, x::Real) = gammalogpdf(promote(k, θ, x)...)
function gammalogpdf(k::T, θ::T, x::T) where {T <: Real}
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

function gammacdf(k::T, θ::T, x::T) where {T <: Real}
    # Handle the degenerate case
    if iszero(k)
        return float(last(promote(x, x >= 0)))
    end
    return first(gamma_inc(k, max(0, x) / θ))
end
gammacdf(k::Real, θ::Real, x::Real) = gammacdf(map(float, promote(k, θ, x))...)

function gammaccdf(k::T, θ::T, x::T) where {T <: Real}
    # Handle the degenerate case
    if iszero(k)
        return float(last(promote(x, x < 0)))
    end
    return last(gamma_inc(k, max(0, x) / θ))
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

    xdθ = max(0, x) / θ
    l, u = gamma_inc(k, xdθ)
    if l < floatmin(Float64) && isfinite(k) && isfinite(xdθ)
        return -log(k) + k * log(xdθ) - xdθ + log(_₁F₁(1.0, 1 + k, xdθ)) - loggamma(k)
    elseif l < 0.7
        return log(l)
    else
        return log1p(-u)
    end
end
function _gammalogcdf(k::T, θ::T, x::T) where {T <: Union{Float16, Float32}}
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

    xdθ = max(0, x) / θ
    l, u = gamma_inc(k, xdθ)
    if u < floatmin(Float64)
        return loggamma(k, xdθ) - loggamma(k)
    elseif u < 0.7
        return log(u)
    else
        return log1p(-l)
    end
end

function _gammalogccdf(k::T, θ::T, x::T) where {T <: Union{Float16, Float32}}
    return T(_gammalogccdf(Float64(k), Float64(θ), Float64(x)))
end

function gammainvcdf(k::Real, θ::Real, p::Real)
    _k, _θ, _p = promote(k, θ, p)
    return _θ * gamma_inc_inv(_k, _p, 1 - _p)
end

function gammainvccdf(k::Real, θ::Real, p::Real)
    _k, _θ, _p = promote(k, θ, p)
    return _θ * gamma_inc_inv(_k, 1 - _p, _p)
end

# Newton iteration in log-space for inverting log-CDF.
# Solves gammalogcdf(k, θ, x) = lp for x.
function _gammainvlogcdf(k::Float64, θ::Float64, lp::Float64)
    if lp > -1
        # Use -expm1(lp) = 1-exp(lp) for accuracy near lp=0
        return Float64(gammainvccdf(k, θ, -expm1(lp)))
    end

    # Use Newton iteration in log-space.
    # Initial estimate from small-x asymptotic:
    # P(k, x/θ) ≈ (x/θ)^k / (k·Γ(k)) for small x/θ
    # ⟹ x₀/θ ≈ exp((lp + log(k) + loggamma(k)) / k)
    logz = (lp + log(k) + loggamma(k)) / k
    if logz + log(θ) < log(nextfloat(0.0))
        return 0.0  # answer underflows Float64
    end
    z = exp(logz)
    x = z * θ
    x = max(x, nextfloat(0.0))

    for _ in 1:200
        lcdf = Float64(gammalogcdf(k, θ, x))
        residual = lcdf - lp
        if abs(residual) <= 2 * eps(max(abs(lp), 1.0))
            break
        end
        lpdf = Float64(gammalogpdf(k, θ, x))
        # Newton step: Δx = residual * F(x)/f(x) = residual * exp(logcdf - logpdf)
        dx = residual * exp(lcdf - lpdf)
        # Guard: don't let x go negative or jump too far
        x_new = x - clamp(dx, -x / 2, x / 2)
        x_new = max(x_new, nextfloat(0.0))
        if x_new == x
            break
        end
        x = x_new
    end

    return x
end

# Newton iteration in log-space for inverting log-CCDF.
# Solves gammalogccdf(k, θ, x) = lp for x.
function _gammainvlogccdf(k::Float64, θ::Float64, lp::Float64)
    if lp > -1
        # Use -expm1(lp) = 1-exp(lp) for accuracy near lp=0
        return Float64(gammainvcdf(k, θ, -expm1(lp)))
    end

    # Extreme right tail: initial estimate
    # Q(k, x/θ) ≈ (x/θ)^(k-1)·exp(-x/θ)/Γ(k) for large x/θ
    # ⟹ log(Q) ≈ (k-1)·log(x/θ) - x/θ - loggamma(k)
    # For large x/θ, the -x/θ term dominates: x₀/θ ≈ -lp - loggamma(k)
    z = max(-lp - loggamma(k), 1.0)
    x = z * θ

    for _ in 1:200
        lccdf = Float64(gammalogccdf(k, θ, x))
        residual = lccdf - lp
        if abs(residual) <= 2 * eps(max(abs(lp), 1.0))
            break
        end
        lpdf = Float64(gammalogpdf(k, θ, x))
        # d/dx log(ccdf(x)) = -f(x)/ccdf(x) = -exp(logpdf - logccdf)
        # Newton step: Δx = (logccdf - lp) / (-exp(logpdf - logccdf))
        #            = -(logccdf - lp) * exp(logccdf - logpdf)
        dx = -residual * exp(lccdf - lpdf)
        # Guard: don't overshoot
        x_new = x - clamp(dx, -x, x / 2)
        x_new = max(x_new, nextfloat(0.0))
        if x_new == x
            break
        end
        x = x_new
    end

    return x
end

function gammainvlogcdf(k::Real, θ::Real, lp::Real)
    T = float(Base.promote_typeof(k, θ, lp))
    _lp = Float64(lp)

    if isnan(_lp) || _lp > 0
        return convert(T, NaN)
    elseif _lp == 0
        return convert(T, Inf)
    elseif isinf(_lp)  # -Inf
        return zero(T)
    end

    _k = Float64(k); _θ = Float64(θ)
    if iszero(_k)
        return _lp >= 0 ? zero(T) : convert(T, NaN)
    end

    return convert(T, _gammainvlogcdf(_k, _θ, _lp))
end

function gammainvlogccdf(k::Real, θ::Real, lp::Real)
    T = float(Base.promote_typeof(k, θ, lp))
    _lp = Float64(lp)

    if isnan(_lp) || _lp > 0
        return convert(T, NaN)
    elseif _lp == 0
        return zero(T)
    elseif isinf(_lp)  # -Inf
        return convert(T, Inf)
    end

    _k = Float64(k); _θ = Float64(θ)
    if iszero(_k)
        return _lp >= 0 ? zero(T) : convert(T, NaN)
    end

    return convert(T, _gammainvlogccdf(_k, _θ, _lp))
end
