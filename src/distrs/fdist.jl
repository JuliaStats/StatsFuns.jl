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
fdistlogpdf(ν1::Real, ν2::Real, x::Number) = (ν1 * log(ν1 * x) + ν2 * log(ν2) - (ν1 + ν2) * log(ν1 * x + ν2)) / 2 - log(x) - logbeta(ν1 / 2, ν2 / 2)

# ChainRules adjoints
ChainRulesCore.@scalar_rule(
    fdistlogpdf(ν1::Real, ν2::Real, x::Number),
    @setup(
        xν1 = x * ν1,
        temp1 = xν1 + ν2,
        a = (x - 1) / temp1,
        νsum = ν1 + ν2,
        di = digamma(νsum / 2),
    ),
    (
        (-log1p(ν2 / xν1) - ν2 * a + di - digamma(ν1 / 2)) / 2,
        (-log1p(xν1 / ν2) + ν1 * a + di - digamma(ν2 / 2)) / 2,
        ((ν1 - 2) / x - ν1 * νsum / temp1) / 2,
    ),
)
