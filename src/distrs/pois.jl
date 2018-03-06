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

function poislogpdf(λ::Real, x::Real)
    (λ ≥ 0 && x ≥ 0) || throw(ArgumentError("λ = $λ and x = $x must both be ≥ 0"))
    T = typeof(float(λ) * x)
    λ = T(λ)
    x = T(x)
    isinteger(x) || return T(-Inf)
    iszero(λ) && return iszero(x) ? zero(T) : T(-Inf)
    iszero(x) && return -λ
    x * log(λ) - λ - lgamma(x + one(x))
end
