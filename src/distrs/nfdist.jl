# functions related to noncentral F distribution

# R implementations
using .RFunctions:
    nfdistpdf,
    nfdistlogpdf,
    nfdistcdf,
    nfdistccdf,
    nfdistlogcdf,
    nfdistlogccdf,
    nfdistinvcdf,
    nfdistinvccdf,
    nfdistinvlogcdf,
    nfdistinvlogccdf

nfdistlogupdf(k1::Real, k2::Real, λ::Real, x::Real) = nfdistlogpdf(k1, k2, λ, x)
nfdistlogulikelihood(k1::Real, k2::Real, λ::Real, x::Real) = nfdistlogpdf(k1, k2, λ, x)
