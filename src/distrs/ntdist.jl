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
    ONE = one(T)
    Av = _₁F₁((k + 1) / 2, ONE / 2, (λ * x)^2 / (2 * muladd(x, x, k)))
    Bv = sqrt(2*ONE) * λ * x * gamma(k/2 +1) / (sqrt(muladd(x, x, k)) * gamma((k+1)/2)) * 
        _₁F₁(k/2 + 1, 3*ONE/2, (λ * x)^2 / (2 * muladd(x, x, k)))
    tdistlogpdf(k, x) - λ^2 / 2 +
    log(Av + Bv)
end
