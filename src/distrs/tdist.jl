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
function tdistpdf(ν::T, x::T) where {T<:Real}
    gamma((ν + 1) / 2) / (sqrt(ν * pi) * gamma(ν / 2)) * (1 + x^2 / ν)^(-(ν + 1) / 2)
end

tdistpdf(ν::Real, x::Real) = tdistpdf(promote(ν, x)...)

# logpdf for numbers with generic types
function tdistlogpdf(ν::T, x::T) where {T<:Real}
    lgamma((ν + 1) / 2) - log(ν * pi) / 2 - lgamma(ν / 2) + (-(ν + 1) / 2) * log(1 + x^2 / ν)
end

tdistlogpdf(ν::Real, x::Real) = tdistlogpdf(promote(ν, x)...)
