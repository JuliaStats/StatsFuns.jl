# functions related to Laplace distribution

laplacepdf(z::Number) = exp(-abs(z)) / 2
laplacepdf(μ::Real, θ::Real, x::Number) = laplacepdf(zval(μ, θ, x)) / θ

laplacelogpdf(z::Number) = -(abs(z) + logtwo)
# Promote to ensure that we don't compute `log(θ)` in low precision and then promote
laplacelogpdf(μ::Real, θ::Real, x::Real) = laplacelogpdf(promote(μ, θ, x)...)
laplacelogpdf(μ::T, θ::T, x::T) where {T <: Real} = laplacelogpdf(zval(μ, θ, x)) - log(θ)

laplacecdf(z::Number) = z < 0 ? exp(z) / 2 : 1 - exp(-z) / 2
laplacecdf(μ::Real, θ::Real, x::Number) = laplacecdf(zval(μ, θ, x))

laplaceccdf(z::Number) = laplacecdf(-z)
laplaceccdf(μ::Real, θ::Real, x::Number) = laplaceccdf(zval(μ, θ, x))

laplacelogcdf(z::Number) = z < 0 ? z - logtwo : log1p(-exp(-z) / 2)
laplacelogcdf(μ::Real, θ::Real, x::Number) = laplacelogcdf(zval(μ, θ, x))

laplacelogccdf(z::Number) = laplacelogcdf(-z)
laplacelogccdf(μ::Real, θ::Real, x::Number) = laplacelogccdf(zval(μ, θ, x))

laplaceinvcdf(p::Real) = p < 1 / 2 ? log(2 * p) : -log(2 * (1 - p))
# Promote to ensure that we don't compute the quantile in low precision and then promote
laplaceinvcdf(μ::Real, θ::Real, p::Real) = laplaceinvcdf(promote(μ, θ, p)...)
laplaceinvcdf(μ::T, θ::T, p::T) where {T <: Real} = xval(μ, θ, laplaceinvcdf(p))

laplaceinvccdf(p::Real) = -laplaceinvcdf(p)
laplaceinvccdf(μ::Real, θ::Real, p::Real) = laplaceinvccdf(promote(μ, θ, p)...)
laplaceinvccdf(μ::T, θ::T, p::T) where {T <: Real} = xval(μ, θ, laplaceinvccdf(p))

function laplaceinvlogcdf(lp::Real)
    x = float(lp)
    # `log(2p)`; the low part of `log(2)` makes its sign reliable at the median
    y = (x + oftype(x, logtwo)) + _logtwo_lo(typeof(x))
    return if y <= 0
        y
    elseif -2 * x >= logtwo
        -log1p(-expm1(y))
    else
        -log(-2 * expm1(x))
    end
end
laplaceinvlogcdf(μ::Real, θ::Real, lp::Real) = laplaceinvlogcdf(promote(μ, θ, lp)...)
laplaceinvlogcdf(μ::T, θ::T, lp::T) where {T <: Real} = xval(μ, θ, laplaceinvlogcdf(lp))

laplaceinvlogccdf(lp::Real) = -laplaceinvlogcdf(lp)
laplaceinvlogccdf(μ::Real, θ::Real, lp::Real) = laplaceinvlogccdf(promote(μ, θ, lp)...)
laplaceinvlogccdf(μ::T, θ::T, lp::T) where {T <: Real} = xval(μ, θ, laplaceinvlogccdf(lp))

# `log(2) = T(logtwo) + _logtwo_lo(T)`
_logtwo_lo(::Type{T}) where {T <: Real} = T(BigFloat(logtwo; precision = 2 * precision(T) + 32) - T(logtwo))
for T in (Float16, Float32, Float64)
    @eval _logtwo_lo(::Type{$T}) = $(_logtwo_lo(T))
end
