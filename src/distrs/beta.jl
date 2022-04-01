# functions related to beta distributions

using HypergeometricFunctions: _₂F₁

# R implementations
using .RFunctions:
    # betapdf,
    # betalogpdf,
    # betacdf,
    # betaccdf,
    # betalogcdf,
    # betalogccdf,
    # betainvcdf,
    # betainvccdf,
    betainvlogcdf,
    betainvlogccdf

# Julia implementations
betapdf(α::Real, β::Real, x::Real) = exp(betalogpdf(α, β, x))

betalogpdf(α::Real, β::Real, x::Real) = betalogpdf(promote(α, β, x)...)
function betalogpdf(α::T, β::T, x::T) where {T<:Real}
    # we ensure that `log(x)` and `log1p(-x)` do not error
    y = clamp(x, 0, 1)
    val = xlogy(α - 1, y) + xlog1py(β - 1, -y) - logbeta(α, β)
    return x < 0 || x > 1 ? oftype(val, -Inf) : val
end

function betacdf(α::Real, β::Real, x::Real)
    # Handle a degenerate case
    if iszero(α) && β > 0
        return float(last(promote(α, β, x, x >= 0)))
    end

    return first(beta_inc(α, β, clamp(x, 0, 1)))
end

function betaccdf(α::Real, β::Real, x::Real)
    # Handle a degenerate case
    if iszero(α) && β > 0
        return float(last(promote(α, β, x, x < 0)))
    end

    last(beta_inc(α, β, clamp(x, 0, 1)))
end

# The log version is currently based on non-log version. When the cdf is very small we shift
# to an implementation based on the hypergeometric function ₂F₁ to avoid underflow.
function betalogcdf(α::T, β::T, x::T) where {T<:Real}
    # Handle a degenerate case
    if iszero(α) && β > 0
        return log(last(promote(x, x >= 0)))
    end

    _x = clamp(x, 0, 1)
    p, q = beta_inc(α, β, _x)
    if p < floatmin(p)
        # see https://dlmf.nist.gov/8.17#E7
        return -log(α) + xlogy(α, _x) + log(_₂F₁(promote(α, 1 - β, α + 1, _x)...)) - logbeta(α, β)
    elseif p <= 0.7
        return log(p)
    else
        return log1p(-q)
    end
end
betalogcdf(α::Real, β::Real, x::Real) = betalogcdf(promote(α, β, x)...)

function betalogccdf(α::Real, β::Real, x::Real)
    # Handle a degenerate case
    if iszero(α) && β > 0
        return log(last(promote(α, β, x, x < 0)))
    end

    p, q = beta_inc(α, β, clamp(x, 0, 1))
    if q < 0.7
        return log(q)
    else
        return log1p(-p)
    end
end

betainvcdf(α::Real, β::Real, p::Real) = first(beta_inc_inv(α, β, p))

betainvccdf(α::Real, β::Real, p::Real) = last(beta_inc_inv(β, α, p))
