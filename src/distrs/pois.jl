# functions related to Poisson distribution

# R implementations
using .RFunctions:
    # poispdf,
    # poislogpdf,
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
function poislogpdf(λ::T, x::T) where {T <: Real}
    val = xlogy(x, λ) - λ - loggamma(x + 1)
    return x >= 0 && isinteger(x) ? val : oftype(val, -Inf)
end
