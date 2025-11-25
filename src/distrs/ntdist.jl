# functions related to noncentral T distribution

# R implementations
using .RFunctions:
    ntdistpdf,
    ntdistlogpdf,
    ntdistcdf,
    ntdistccdf,
    ntdistlogcdf,
    ntdistlogccdf,
    ntdistinvcdf,
    ntdistinvccdf,
    ntdistinvlogcdf,
    ntdistinvlogccdf

ntdistlogupdf(k::Real, λ::Real, x::Real) = ntdistlogpdf(k, λ, x)
ntdistlogulikelihood(k::Real, λ::Real, x::Real) = ntdistlogpdf(k, λ, x)
