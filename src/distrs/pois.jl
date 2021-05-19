# functions related to Poisson distribution

import .RFunctions:
    # poispdf,
    # poislogpdf,
    # poiscdf,
    # poisccdf,
    # poislogcdf,
    # poislogccdf,
    poisinvcdf,
    poisinvccdf,
    poisinvlogcdf,
    poisinvlogccdf

poispdf(λ::Real, x::Real) = exp(poislogpdf(λ, x))

poislogpdf(λ::T, x::T) where {T <: Real} = xlogy(x, λ) - λ - loggamma(x + 1)
poislogpdf(λ::Real, x::Real) = poislogpdf(promote(float(λ), x)...)

# Just use the Gamma definitions
poiscdf(λ::Real, x::Real) = gammaccdf(floor(x + 1), 1, λ)

poisccdf(λ::Real, x::Real) = gammacdf(floor(x + 1), 1, λ)

poislogcdf(λ::Real, x::Real) = gammalogccdf(floor(x + 1), 1, λ)

poislogccdf(λ::Real, x::Real) = gammalogcdf(floor(x + 1), 1, λ)
