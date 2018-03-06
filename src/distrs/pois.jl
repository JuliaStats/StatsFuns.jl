# functions related to Poisson distribution

import .RFunctions:
    poispdf,
    poislogpdf,
    poiscdf,
    poisccdf,
    poislogcdf,
    poislogccdf,
    poisinvcdf,
    poisinvccdf,
    poisinvlogcdf,
    poisinvlogccdf

# generic versions
poispdf(λ::Real, x::Real) = exp(poislogpdf(λ, x))

function poislogpdf(λ::T, x::T) where {T <: Real}
    iszero(λ) ? (iszero(x) ? x : T(-Inf)) : x * log(λ) - λ - lgamma(x + 1)
end

poislogpdf(λ::Number, x::Number) = poislogpdf(promote(float(λ), x)...)
