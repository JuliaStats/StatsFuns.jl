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

poislogpdf(λ::T, x::T) where {T <: Real} = xlogy(x, λ) - λ - loggamma(x + 1)

poislogpdf(λ::Number, x::Number) = poislogpdf(promote(float(λ), x)...)

function poislogcdf(λ::T, x::T) where {T<:Real}
    return loggamma(floor(x + 1), float(λ)) - logfactorial(floor(x))
end
poislogcdf(λ::Number, x::Number) = poislogcdf(promote(float(λ), x)...)

poiscdf(λ::T, x::T; precision=0) where {T<:Real} = last(gamma_inc(floor(x + 1), λ, precision))
poiscdf(λ::Number, x::Number; precision=0) = poiscdf(promote(float(λ), x)...)
poisccdf(λ::T, x::T; precision=0) where {T<:Real} = first(gamma_inc(floor(x + 1), λ, precision))
poisccdf(λ::Number, x::Number; precision=0) = poisccdf(promote(float(λ), x)...)

#=
function poislogpdf(λ::Union{Float32,Float64}, x::Union{Float64,Float32,Integer})
    if iszero(λ)
        iszero(x) ? zero(λ) : oftype(λ, -Inf)
    elseif iszero(x)
        -λ
    else
        -lstirling_asym(x + 1)
=#
