# functions related to beta distributions

import .Rmath:
    betalogpdf,
    betacdf,
    betaccdf,
    betalogcdf,
    betalogccdf,
    betainvcdf,
    betainvccdf,
    betainvlogcdf,
    betainvlogccdf

# pdf
betapdf(α::Float64, β::Float64, x::Float64) = exp(betalogpdf(α, β, x))
betapdf(α::Real, β::Real, x::Real) = betapdf(f64(α), f64(β), f64(x))
