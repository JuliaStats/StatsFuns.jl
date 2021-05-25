# functions related to F distribution

import .RFunctions:
    # fdistpdf,
    # fdistlogpdf,
    fdistcdf,
    fdistccdf,
    fdistlogcdf,
    fdistlogccdf,
    fdistinvcdf,
    fdistinvccdf,
    fdistinvlogcdf,
    fdistinvlogccdf

# pdf for numbers with generic types
fdistpdf(ν1::Real, ν2::Real, x::Real) = exp(fdistlogpdf(ν1, ν2, x))

# logpdf for numbers with generic types
function fdistlogpdf(ν1::T, ν2::T, x::T) where T<:Real
    ν1ν2 = ν1/ν2
    return  xlogy(ν1/2, ν1ν2) + xlogy(ν1/2 - 1, x) - xlogy((ν1 + ν2)/2, 1 + ν1ν2*x) - logbeta(ν1/2, ν2/2)
end
fdistlogpdf(ν1::Real, ν2::Real, x::Real) = fdistlogpdf(promote(ν1, ν2, x)...)
