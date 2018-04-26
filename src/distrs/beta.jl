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
function betapdf(α::T, β::T, x::V) where {T <: Real, V <: Number}
    if 0 <= x <= 1
       x^(α - 1) * (1 - x)^(β - 1) / beta(α, β) 
    else
       zero(T)
    end
end
# logpdf for numbers with generic types
function betalogpdf(α::T, β::T, x::V) where {T <: Real, V <: Number}
    if 0 <= x <= 1
        (α - 1) * log(x) + (β - 1) * log1p(-x) - lbeta(α, β)
    else
       convert(T, -Inf)
    end
end
