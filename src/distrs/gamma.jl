# functions related to gamma distribution

import .RFunctions:
    gammapdf,
    gammalogpdf,
    gammacdf,
    gammaccdf,
    gammalogcdf,
    gammalogccdf,
    gammainvcdf,
    gammainvccdf,
    gammainvlogcdf,
    gammainvlogccdf

# pdf for numbers with generic types
gammapdf(k::Real, θ::Real, x::Real) = exp(gammalogpdf(k, θ, x))

# logpdf for numbers with generic types
function gammalogpdf(k::T, θ::T, x::T) where {T<:Real}
    k > 0 && θ > 0 || throw(ArgumentError("shape k = $k and scale θ = $θ must both be positive"))
    x > 0 ? -lgamma(k) - k * log(θ) + (k - 1) * log(x) - x / θ : T(-Inf)
end

gammalogpdf(k::Real, θ::Real, x::Real) = gammalogpdf(promote(k, θ, x)...)
