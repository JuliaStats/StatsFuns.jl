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

function poislogpdf(λ::T, x::T) where {T <: AbstractFloat}
    λ < 0 && throw(ArgumentError("λ = $λ must be non-negative"))
    isinteger(x) && 0 ≤ x || return T(-Inf)
    iszero(λ) && return iszero(x) ? zero(T) : T(-Inf)
    xlogy(x, λ) - λ - lgamma(x + 1)
end    

poislogpdf(λ::Real, x::Real) = poislogpdf(promote(float(λ), x)...)
