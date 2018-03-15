# functions related to beta distributions

import .RFunctions:
    betapdf,
    betalogpdf,
    betacdf,
    betaccdf,
    betalogcdf,
    betalogccdf,
    betainvcdf,
    betainvccdf,
    betainvlogcdf,
    betainvlogccdf

# pdf for numbers with generic types
betapdf(α::Real, β::Real, x::Real) = exp(betalogpdf(α, β, x))

# logpdf for numbers with generic types
function betalogpdf(α::T, β::T, x::T) where {T<:Real}
    α > 0 && β > 0 || throw(ArgumentError("α and β must both be positive, got α = $α, β = $β "))
    0 < x < 1 ? (α - 1) * log(x) + (β - 1) * log1p(-x) - lbeta(α, β) : T(-Inf)
end

betalogpdf(α::Real, β::Real, x::Real) = betalogpdf(promote(float(α), β, x)...)
