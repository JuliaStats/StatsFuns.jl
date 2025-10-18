# functions related to noncentral chi-square distribution

# R implementations
using .RFunctions:
    nchisqpdf,
    nchisqlogpdf,
    nchisqcdf,
    nchisqccdf,
    nchisqlogcdf,
    nchisqlogccdf,
    nchisqinvcdf,
    nchisqinvccdf,
    nchisqinvlogcdf,
    nchisqinvlogccdf

nchisqlogupdf(k::Real, λ::Real, x::Real) = nchisqlogpdf(k, λ, x)
nchisqlogulikelihood(k::Real, λ::Real, x::Real) = nchisqlogpdf(k, λ, x)
