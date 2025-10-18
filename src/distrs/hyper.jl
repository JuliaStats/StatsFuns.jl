# functions related to hyper-geometric distribution

# R implementations
using .RFunctions:
    hyperpdf,
    hyperlogpdf,
    hypercdf,
    hyperccdf,
    hyperlogcdf,
    hyperlogccdf,
    hyperinvcdf,
    hyperinvccdf,
    hyperinvlogcdf,
    hyperinvlogccdf


hyperlogupdf(ms::Real, mf::Real, n::Real, x::Real) = hyperlogpdf(ms, mf, n, x)
hyperlogulikelihood(ms::Real, mf::Real, n::Real, x::Real) = hyperlogpdf(ms, mf, n, x)
