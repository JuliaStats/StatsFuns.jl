# functions related to noncentral T distribution

# R implementations
using .RFunctions:
    ntdistcdf,
    ntdistccdf,
    ntdistlogcdf,
    ntdistlogccdf,
    ntdistinvcdf,
    ntdistinvccdf,
    ntdistinvlogcdf,
    ntdistinvlogccdf

using SpecialFunctions: gamma

function ntdistpdf(k::Real, λ::Real, x::Real)
    return exp(ntdistlogpdf(k, λ, x))
end

ntdistlogpdf(k::Real, λ::Real, x::Real) = ntdistlogpdf(promote(k, λ, x)...)

function ntdistlogpdf(k::T, λ::T, x::T) where {T<:Real}
    tdistlogpdf(k, x) - λ^2 / 2 +
    log(_₁F₁((k + 1) / 2, 1 / 2, (λ * x)^2 / (2 * muladd(x,x,k))) +
    sqrt(2) * λ * x * gamma(k/2 +1) / (sqrt(muladd(x,x,k))*gamma((k+1)/2)) * 
        _₁F₁(k/2 + 1, 3/2, (λ * x)^2 / (2 * muladd(x,x,k)))
    )
end
