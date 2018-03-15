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

poislogpdf(λ::T, x::T) where {T <: Real} = xlogy(x, λ) - λ - lgamma(x + 1)

poislogpdf(λ::Number, x::Number) = poislogpdf(promote(float(λ), x)...)

#=
function poislogpdf(λ::Union{Float32,Float64}, x::Union{Float64,Float32,Integer})
    if iszero(λ)
        iszero(x) ? zero(λ) : oftype(λ, -Inf)
    elseif iszero(x)
        -λ
    else
        -lstirling_asym(x + 1)
=#