# functions related to negative binomial distribution

# R implementations
using .RFunctions:
    nbinompdf,
    nbinomlogpdf,
    nbinomcdf,
    nbinomccdf,
    nbinomlogcdf,
    nbinomlogccdf,
    nbinominvcdf,
    nbinominvccdf,
    nbinominvlogcdf,
    nbinominvlogccdf

nbinomlogupdf(r::Real, p::Real, x::Real) = nbinomlogpdf(r, p, x)
nbinomlogulikelihood(r::Real, p::Real, x::Real) = nbinomlogpdf(r, p, x)
