# functions related to gamma distribution

import .RFunctions:
    gammacdf,
    gammaccdf,
    gammalogcdf,
    gammalogccdf,
    gammainvcdf,
    gammainvccdf,
    gammainvlogcdf,
    gammainvlogccdf


# pdf
gammapdf(k::Real, θ::Real, x::Number) = 1 / (gamma(k) * θ^k) * x^(k - 1) * exp(-x / θ)

# logpdf
gammalogpdf(k::Real, θ::Real, x::Number) = -lgamma(k) - k * log(θ) + (k - 1) * log(x) - x / θ
