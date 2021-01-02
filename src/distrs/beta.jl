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
betapdf(α::Real, β::Real, x::Number) = x^(α - 1) * (1 - x)^(β - 1) / beta(α, β)

# logpdf for numbers with generic types
betalogpdf(α::Real, β::Real, x::Number) = (α - 1) * log(x) + (β - 1) * log1p(-x) - logbeta(α, β)

# ChainRules adjoint
ChainRulesCore.@scalar_rule(
    betalogpdf(α::Real, β::Real, x::Number),
    @setup(z = digamma(α + β)),
    (
        log(x) + z - digamma(α),
        log1p(-x) + z - digamma(β),
        (α - 1) / x + (1 - β) / (1 - x),
    ),
)
