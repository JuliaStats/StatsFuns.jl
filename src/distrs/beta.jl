# functions related to beta distributions

using HypergeometricFunctions: _₂F₁

import .RFunctions:
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

# pdf for numbers with generic types
betapdf(α::Real, β::Real, x::Real) = exp(betalogpdf(α, β, x))

# logpdf for numbers with generic types
function betalogpdf(α::T, β::T, x::T) where T<:Real
    return if x <= 0.5
        xlogy(α - 1, x) + (β - 1) * log1p(-x) - logbeta(α, β)
    else
        betalogpdf(β, α, 1 - x)
    end
end
betalogpdf(α::Real, β::Real, x::Real) = betalogpdf(promote(α, β, x)...)

betacdf(α::Real, β::Real, x::Real) = first(beta_inc(α, β, x))

betaccdf(α::Real, β::Real, x::Real) = last(beta_inc(α, β, x))

# The log version is currently based on non-log version. When the cdf is very small we shift
# to a implementation based on the hypergeometric function ₂F₁ to avoid underflow.
function betalogcdf(α::Float64, β::Float64, x::Float64)
    p, q = beta_inc(α, β, x)
    if p < eps(one(p))
        # see https://dlmf.nist.gov/8.17#E7
        return -log(α) + α*log(x) + log(_₂F₁(promote(α, 1 - β, α + 1, x)...)) - logbeta(α, β)
    elseif p <= 0.7
        return log(p)
    else
        return log1p(-q)
    end
end
betalogcdf(α::Real, β::Real, x::Real) = betalogcdf(promote(float(α), β, x)...)
betalogcdf(k::T, x::T) where T = throw(MethodError(betalogcdf, (α, β, x)))

function betalogccdf(α::Real, β::Real, x::Real)
    p, q = beta_inc(α, β, x)
    if q < 0.7
        return log(q)
    else
        return log1p(-p)
    end
end

betainvcdf(α::Float64, β::Float64, p::Float64) = first(beta_inc_inv(α, β, p, 1 - p))
betainvcdf(α::Real, β::Real, p::Real) = betainvcdf(promote(float(α), β, p)...)
betainvcdf(α::T, β::T, p::T) where T = throw(MethodError(betainvcdf, (α, β, p)))

betainvccdf(α::Float64, β::Float64, p::Float64) = first(beta_inc_inv(α, β, 1 - p, p))
betainvccdf(α::Real, β::Real, p::Real) = betainvccdf(promote(float(α), β, p)...)
betainvccdf(α::T, β::T, p::T) where T = throw(MethodError(betainvccdf, (α, β, p)))
