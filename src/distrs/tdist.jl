# functions related to student's T distribution

# R implementations
using .RFunctions:
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

# Julia implementations
tdistpdf(ν::Real, x::Real) = exp(tdistlogpdf(ν, x))

tdistlogpdf(ν::Real, x::Real) = tdistlogpdf(promote(ν, x)...)
function tdistlogpdf(ν::T, x::T) where {T<:Real}
    isinf(ν) && return normlogpdf(x)
    νp12 = (ν + 1) / 2
    return loggamma(νp12) - (logπ + log(ν)) / 2 - loggamma(ν / 2) - νp12 * log1p(x^2 / ν)
end
