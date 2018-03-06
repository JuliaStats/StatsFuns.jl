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
    λ = float(λ)
    nInf = oftype(λ, -Inf) * one(x)
    isinteger(x) || return nInf
    iszero(λ) && return (iszero(x) ? zero(λ) : nInf) * one(x)
    iszero(x) && return -λ * one(x)
    x * log(λ) - λ - lgamma(one(λ) * (x + one(x)))
end

# The peculiar argument in `lgamma(one(λ) * (x + one(x)))` is to retain precision
# when `λ` is a BigFloat but also allow for `x` to be a Dual number.
