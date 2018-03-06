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
function fdistpdf(d1::T, d2::T, x::T) where {T<:Real}
    sqrt((d1 * x)^d1 * d2^d2 / (d1 * x + d2)^(d1 + d2)) / (x * beta(d1 / 2, d2 / 2))
end

fdistpdf(d1::Real, d2::Real, x::Real) = fdistpdf(promote(d1, d2, x)...)

# logpdf for numbers with generic types
function fdistlogpdf(d1::T, d2::T, x::T) where {T<:Real}
    (d1 * log(d1 * x) + d2 * log(d2) - (d1 + d2) * log(d1 * x + d2)) / 2 - log(x) - lbeta(d1 / 2, d2 / 2)
end

fdistlogpdf(d1::Real, d2::Real, x::Real) = fdistlogpdf(promote(d1, d2, x)...) 