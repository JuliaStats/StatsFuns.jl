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
tdistpdf(ν::Real, x::Real) = exp(tdistlogpdf(ν, x))

# logpdf for numbers with generic types
function tdistlogpdf(ν::T, x::T) where {T<:Real}
    ν > 0 || throw(ArgumentError("ν should be positive, got ν = $ν"))
    lgamma((ν + 1) / 2) - log(ν * pi) / 2 - lgamma(ν / 2) + (-(ν + 1) / 2) * log(1 + abs2(x) / ν)
end

tdistlogpdf(ν::Real, x::Real) = tdistlogpdf(promote(ν, float(x))...)
