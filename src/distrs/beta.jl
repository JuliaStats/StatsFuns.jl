# functions related to beta distributions

using HypergeometricFunctions: _₂F₁

# Julia implementations
betapdf(α::Real, β::Real, x::Real) = exp(betalogpdf(α, β, x))

betalogpdf(α::Real, β::Real, x::Real) = betalogpdf(promote(α, β, x)...)
function betalogpdf(α::T, β::T, x::T) where {T <: Real}
    # we ensure that `log(x)` and `log1p(-x)` do not error
    y = clamp(x, 0, 1)
    val = xlogy(α - 1, y) + xlog1py(β - 1, -y) - logbeta(α, β)
    return x < 0 || x > 1 ? oftype(val, -Inf) : val
end

function betacdf(α::Real, β::Real, x::Real)
    # Handle degenerate cases
    if iszero(α) && β > 0
        return float(last(promote(α, β, x, x >= 0)))
    elseif iszero(β) && α > 0
        return float(last(promote(α, β, x, x >= 1)))
    end

    return first(beta_inc(α, β, clamp(x, 0, 1)))
end

function betaccdf(α::Real, β::Real, x::Real)
    # Handle degenerate cases
    if iszero(α) && β > 0
        return float(last(promote(α, β, x, x < 0)))
    elseif iszero(β) && α > 0
        return float(last(promote(α, β, x, x < 1)))
    end

    return last(beta_inc(α, β, clamp(x, 0, 1)))
end

# The log version is currently based on non-log version. When the cdf is very small we shift
# to an implementation based on the hypergeometric function ₂F₁ to avoid underflow.
function betalogcdf(α::T, β::T, x::T) where {T <: Real}
    # Handle degenerate cases
    if iszero(α) && β > 0
        return log(last(promote(x, x >= 0)))
    elseif iszero(β) && α > 0
        return log(last(promote(x, x >= 1)))
    end

    _x = clamp(x, 0, 1)
    p, q = beta_inc(α, β, _x)
    if p < floatmin(p)
        # see https://dlmf.nist.gov/8.17#E8
        # we use E8 instead of E7 due to https://github.com/JuliaMath/HypergeometricFunctions.jl/issues/47
        return -log(α) + xlogy(α, _x) + xlog1py(β, -_x) + log(_₂F₁(promote(α + β, 1, α + 1, _x)...; method = :positive)) - logbeta(α, β)
    elseif p <= 0.7
        return log(p)
    else
        return log1p(-q)
    end
end
betalogcdf(α::Real, β::Real, x::Real) = betalogcdf(promote(α, β, x)...)

function betalogccdf(α::Real, β::Real, x::Real)
    # Handle degenerate cases
    if iszero(α) && β > 0
        return log(last(promote(α, β, x, x < 0)))
    elseif iszero(β) && α > 0
        return log(last(promote(α, β, x, x < 1)))
    end

    p, q = beta_inc(α, β, clamp(x, 0, 1))
    if q < 0.7
        return log(q)
    else
        return log1p(-p)
    end
end

function betainvcdf(α::Real, β::Real, p::Real)
    # Handle degenerate cases
    if 0 ≤ p ≤ 1
        if iszero(α) && β > 0
            return last(promote(α, β, p, false))
        elseif iszero(β) && α > 0
            return last(promote(α, β, p, p > 0))
        end
    end

    return first(beta_inc_inv(α, β, p))
end

function betainvccdf(α::Real, β::Real, p::Real)
    # Handle degenerate cases
    if 0 ≤ p ≤ 1
        if iszero(α) && β > 0
            return last(promote(α, β, p, p == 0))
        elseif iszero(β) && α > 0
            return last(promote(α, β, p, true))
        end
    end

    return last(beta_inc_inv(β, α, p))
end

# Newton iteration in log-space for inverting log-CDF.
# Solves betalogcdf(α, β, x) = lp for x.
# Uses the identity: d/dx log(F(x)) = f(x)/F(x) = exp(logpdf - logcdf)
# Newton step: x_{n+1} = x_n - (logcdf(x_n) - lp) / exp(logpdf(x_n) - logcdf(x_n))
#            = x_n - (logcdf(x_n) - lp) * exp(logcdf(x_n) - logpdf(x_n))
function _betainvlogcdf(α::Float64, β::Float64, lp::Float64)
    if lp > -1
        # Moderate lp: use -expm1(lp) = 1-exp(lp) for accuracy near lp=0
        return Float64(betainvccdf(α, β, -expm1(lp)))
    end

    # Use Newton iteration in log-space.
    # Initial estimate from small-x asymptotic: I_x(α,β) ≈ x^α / (α·B(α,β))
    # ⟹ log(I_x) ≈ α·log(x) - log(α) - logbeta(α,β)
    # ⟹ x₀ ≈ exp((lp + log(α) + logbeta(α,β)) / α)
    logx = (lp + log(α) + logbeta(α, β)) / α
    if logx < log(nextfloat(0.0))
        return 0.0  # answer underflows Float64
    end
    x = exp(logx)
    x = min(x, 1.0 - eps(1.0))

    for _ in 1:200
        lcdf = Float64(betalogcdf(α, β, x))
        residual = lcdf - lp
        if abs(residual) <= 2 * eps(max(abs(lp), 1.0))
            break
        end
        lpdf = Float64(betalogpdf(α, β, x))
        # Newton step: Δx = residual * F(x)/f(x) = residual * exp(logcdf - logpdf)
        dx = residual * exp(lcdf - lpdf)
        # Guard: don't let x go negative or jump too far
        x_new = x - clamp(dx, -x / 2, x / 2)
        x_new = clamp(x_new, nextfloat(0.0), 1.0 - eps(1.0))
        if x_new == x
            break
        end
        x = x_new
    end

    return x
end

function betainvlogcdf(α::Real, β::Real, lp::Real)
    T = float(Base.promote_typeof(α, β, lp))
    _lp = Float64(lp)

    if isnan(_lp) || _lp > 0
        return convert(T, NaN)
    elseif _lp == 0
        return one(T)
    elseif isinf(_lp)  # -Inf
        return zero(T)
    end

    # Handle degenerate cases
    _α = Float64(α); _β = Float64(β)
    if iszero(_α) && _β > 0
        return zero(T)
    elseif iszero(_β) && _α > 0
        return _lp > log1p(-1.0) ? one(T) : zero(T)
    end

    return convert(T, _betainvlogcdf(_α, _β, _lp))
end

function betainvlogccdf(α::Real, β::Real, lp::Real)
    T = float(Base.promote_typeof(α, β, lp))
    _lp = Float64(lp)

    if isnan(_lp) || _lp > 0
        return convert(T, NaN)
    elseif _lp == 0
        return zero(T)
    elseif isinf(_lp)  # -Inf
        return one(T)
    end

    # Handle degenerate cases
    _α = Float64(α); _β = Float64(β)
    if iszero(_α) && _β > 0
        return _lp > log1p(-1.0) ? zero(T) : one(T)
    elseif iszero(_β) && _α > 0
        return one(T)
    end

    # For moderate lp (near 0): use -expm1(lp) = 1-exp(lp) for accuracy
    if _lp > -1
        return convert(T, Float64(betainvcdf(_α, _β, -expm1(_lp))))
    end

    # For very negative lp: betainvlogccdf(α, β, lp) = 1 - betainvlogcdf(β, α, lp)
    return convert(T, 1.0 - _betainvlogcdf(_β, _α, _lp))
end
