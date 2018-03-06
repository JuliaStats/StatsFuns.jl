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
function gammapdf(k::T, θ::T, x::T) where {T<:Real}
    inv(gamma(k) * θ^k) * x^(k - 1) * exp(-x / θ)
end

gammapdf(k::Real, θ::Real, x::Real) = gammapdf(promote(k, θ, x)...)

# logpdf for numbers with generic types
function gammalogpdf(k::T, θ::T, x::T) where {T<:Real}
    -lgamma(k) - k * log(θ) + (k - 1) * log(x) - x / θ
end

gammalogpdf(k::Real, θ::Real, x::Real) = gammalogpdf(promote(k, θ, x)...)
