# functions related to beta distributions

import .RFunctions:
    betacdf,
    betaccdf,
    betalogcdf,
    betalogccdf,
    betainvcdf,
    betainvccdf,
    betainvlogcdf,
    betainvlogccdf

# pdf
betapdf(α::Real, β::Real, x::Number) = x^(α - 1) * (1 - x)^(β - 1) / beta(α, β)

# logpdf
betalogpdf(α::Real, β::Real, x::Number) = (α - 1) * log(x) + (β - 1) * log(1 - x) - lbeta(α, β)
