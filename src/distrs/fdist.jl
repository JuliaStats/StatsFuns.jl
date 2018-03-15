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
fdistpdf(d1::Real, d2::Real, x::Real) = exp(fdistlogpdf(d1, d2, x))

# logpdf for numbers with generic types
function fdistlogpdf(d1::T, d2::T, x::T) where {T<:Real}
    isinteger(d1) && d1 > 0 && isinteger(d2) && d2 > 0 || 
        throw(ArgumentError("d1 and d2 must be positive integers, got d1 = d2 and d2 = d2"))
    (d1 * log(d1 * x) + d2 * log(d2) - (d1 + d2) * log(d1 * x + d2)) / 2 - log(x) - lbeta(d1 / 2, d2 / 2)
end

fdistlogpdf(d1::Real, d2::Real, x::Real) = fdistlogpdf(promote(d1, d2, float(x))...) 