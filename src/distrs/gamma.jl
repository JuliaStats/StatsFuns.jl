# functions related to gamma distribution

import .RFunctions:
    gammapdf,
    gammalogpdf,
    gammacdf,
    gammaccdf,
    gammalogcdf,
    gammalogccdf,
    gammainvcdf,
    gammainvccdf,
    gammainvlogcdf,
    gammainvlogccdf

# pdf for numbers with generic types
gammapdf(k::Real, θ::Real, x::Number) = 1 / (gamma(k) * θ^k) * x^(k - 1) * exp(-x / θ)

# logpdf for numbers with generic types
gammalogpdf(k::Real, θ::Real, x::Number) = -loggamma(k) - k * log(θ) + (k - 1) * log(x) - x / θ

# ChainRules adjoints
ChainRulesCore.@scalar_rule(
    gammalogpdf(k::Real, θ::Real, x::Number),
    @setup(
        invθ = inv(θ),
        xoθ = invθ * x,
        z = xoθ - k,
    ),
    (
        log(xoθ) - digamma(k),
        invθ * z,
        - (1 + z) / x,
    ),
)
