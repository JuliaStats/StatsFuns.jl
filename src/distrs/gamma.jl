# functions related to gamma distribution

# R implementations
using .RFunctions:
    # gammapdf,
    # gammalogpdf,
    gammacdf,
    gammaccdf,
    gammalogcdf,
    gammalogccdf,
    gammainvcdf,
    gammainvccdf,
    gammainvlogcdf,
    gammainvlogccdf

# Julia implementations
gammapdf(k::Real, θ::Real, x::Real) = exp(gammalogpdf(k, θ, x))

gammalogpdf(k::Real, θ::Real, x::Real) = gammalogpdf(promote(k, θ, x)...)
function gammalogpdf(k::T, θ::T, x::T) where {T<:Real}
    # we ensure that `log(x)` does not error if `x < 0`
    xθ = max(x, 0) / θ
    val = -loggamma(k) + xlogy(k - 1, xθ) - log(θ) - xθ
    return x < 0 ? oftype(val, -Inf) : val
end
