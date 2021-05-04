# functions related to student's T distribution

import .RFunctions:
    # tdistpdf,
    # tdistlogpdf,
    tdistcdf,
    tdistccdf,
    tdistlogcdf,
    tdistlogccdf,
    tdistinvcdf,
    tdistinvccdf,
    tdistinvlogcdf,
    tdistinvlogccdf

# pdf for numbers with generic types
tdistpdf(ν::T, x::T) where T<:Real = gamma((ν + 1) / 2) / (sqrt(ν * π) * gamma(ν / 2)) * (1 + x^2 / ν)^(-(ν + 1) / 2)
tdistpdf(ν::Real, x::Real) = tdistpdf(promote(ν, x)...)

# logpdf for numbers with generic types
tdistlogpdf(ν::T, x::T) where T<:Real = loggamma((ν + 1) / 2) - log(ν * π) / 2 - loggamma(ν / 2) + (-(ν + 1) / 2) * log(1 + x^2 / ν)
tdistlogpdf(ν::Real, x::Real) = tdistlogpdf(promote(ν, x)...)
