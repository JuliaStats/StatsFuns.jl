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
tdistpdf(ν::Real, x::Number) = gamma((ν + 1) / 2) / (sqrt(ν * pi) * gamma(ν / 2)) * (1 + x^2 / ν)^(-(ν + 1) / 2)

# logpdf for numbers with generic types
tdistlogpdf(ν::Real, x::Number) = lgamma((ν + 1) / 2) - log(ν * pi) / 2 - lgamma(ν / 2) + (-(ν + 1) / 2) * log(1 + x^2 / ν)
