# functions related to student's T distribution

import .RFunctions:
    tdistpdf,
    tdistlogpdf,
    tdistcdf,
    tdistccdf,
    tdistlogcdf,
    tdistlogccdf,
    tdistinvcdf,
    tdistinvccdf,
    tdistinvlogcdf,
    tdistinvlogccdf

# pdf for numbers with generic types
tdistpdf(ν::Real, x::Number) = gamma(0.5 * (ν + 1)) / (sqrt(ν * pi) * gamma(0.5 * ν)) * (1 + x^2 / ν)^(-0.5 * (ν + 1))

# logpdf for numbers with generic types
tdistlogpdf(ν::Real, x::Number) = lgamma(0.5 * (ν + 1)) - 0.5 * log(ν * pi) - lgamma(0.5 * ν) + (-0.5 * (ν + 1)) * log(1 + x^2 / ν)
