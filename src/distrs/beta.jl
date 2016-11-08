# functions related to beta distributions

import .RFunctions:
    betapdf,
    betalogpdf,
    betacdf,
    betaccdf,
    betalogcdf,
    betalogccdf,
    betainvcdf,
    betainvccdf,
    betainvlogcdf,
    betainvlogccdf

# pdf for numbers with generic types
betapdf(α::Real, β::Real, x::Real) = x^(α - 1) * (1 - x)^(β - 1) / beta(α, β)

# logpdf for numbers with generic types
betalogpdf(α::Real, β::Real, x::Real) = (α - 1) * log(x) + (β - 1) * log1p(-x) - lbeta(α, β)
