# functions related to Poisson distribution

# R implementations
# For pdf and logpdf we use the Julia implementation
using .RFunctions:
    poiscdf,
    poisccdf,
    poislogcdf,
    poislogccdf,
    poisinvcdf,
    poisinvccdf,
    poisinvlogcdf,
    poisinvlogccdf

# Julia implementations
poispdf(λ::Real, x::Real) = exp(poislogpdf(λ, x))

poislogpdf(λ::Real, x::Real) = poislogpdf(promote(λ, x)...)
poislogpdf(λ::T, x::T) where {T <: Real} = xlogy(x, λ) - λ - loggamma(x + 1)
