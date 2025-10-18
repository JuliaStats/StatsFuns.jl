# functions related to noncentral beta distribution

# R implementations
using .RFunctions:
    nbetapdf,
    nbetalogpdf,
    nbetacdf,
    nbetaccdf,
    nbetalogcdf,
    nbetalogccdf,
    nbetainvcdf,
    nbetainvccdf,
    nbetainvlogcdf,
    nbetainvlogccdf

nbetalogupdf(α::Real, β::Real, λ::Real, x::Real) = nbetalogpdf(α, β, λ, x)
nbetalogulikelihood(α::Real, β::Real, λ::Real, x::Real) = nbetalogpdf(α, β, λ, x)
