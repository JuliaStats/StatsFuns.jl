# functions related to F distribution

import .RFunctions:
    fdistpdf,
    fdistlogpdf,
    fdistcdf,
    fdistccdf,
    fdistlogcdf,
    fdistlogccdf,
    fdistinvcdf,
    fdistinvccdf,
    fdistinvlogcdf,
    fdistinvlogccdf

# pdf for numbers with generic types
fdistpdf(ν1::Real, ν2::Real, x::Number) = sqrt((ν1 * x)^ν1 * ν2^ν2 / (ν1 * x + ν2)^(ν1 + ν2)) / (x * beta(ν1 / 2, ν2 / 2))

# logpdf for numbers with generic types
fdistlogpdf(ν1::Real, ν2::Real, x::Number) = (ν1 * log(ν1 * x) + ν2 * log(ν2) - (ν1 + ν2) * log(ν1 * x + ν2)) / 2 - log(x) - lbeta(ν1 / 2, ν2 / 2)
