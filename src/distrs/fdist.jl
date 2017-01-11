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
fdistpdf(d1::Real, d2::Real, x::Number) = sqrt((d1 * x)^d1 * d2^d2 / (d1 * x + d2)^(d1 + d2)) / (x * beta(d1 / 2, d2 / 2))

# logpdf for numbers with generic types
fdistlogpdf(d1::Real, d2::Real, x::Number) = (d1 * log(d1 * x) + d2 * log(d2) - (d1 + d2) * log(d1 * x + d2)) / 2 - log(x) - lbeta(d1 / 2, d2 / 2)
