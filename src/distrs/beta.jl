# functions related to beta distributions

import .RFunctions:
    # betapdf,
    # betalogpdf,
    betacdf,
    betaccdf,
    betalogcdf,
    betalogccdf,
    betainvcdf,
    betainvccdf,
    betainvlogcdf,
    betainvlogccdf

# pdf for numbers with generic types
betapdf(α::T, β::T, x::T) where T<:Real = x^(α - 1) * (1 - x)^(β - 1) / beta(α, β)
betapdf(α::Real, β::Real, x::Real) = betapdf(promote(α, β, x)...)

# logpdf for numbers with generic types
function betalogpdf(α::T, β::T, x::T) where T<:Real
    return if x <= 0.5
        xlogy(α - 1, x) + (β - 1) * log1p(-x) - logbeta(α, β)
    else
        betalogpdf(β, α, 1 - x)
    end
end
betalogpdf(α::Real, β::Real, x::Real) = betalogpdf(promote(α, β, x)...)
