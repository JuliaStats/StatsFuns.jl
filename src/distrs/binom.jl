# functions related to binomial distribution

import .RFunctions:
    binompdf,
    binomlogpdf,
    binomcdf,
    binomccdf,
    binomlogcdf,
    binomlogccdf,
    binominvcdf,
    binominvccdf,
    binominvlogcdf,
    binominvlogccdf

# pdf for numbers with generic types
binompdf(n::Real, p::Real, k::Real) = exp(binomlogpdf(n, p, k))

# logpdf for numbers with generic types
function binomlogpdf(n::T, p::T, k::T) where {T<:Real}
    isinteger(n) && 0 ≤ n && 0 ≤ p ≤ 1 || 
        throw(ArgumentError("n must be a non-negative integer and 0 ≤ p ≤ 1"))
    isinteger(k) && 0 ≤ k ≤ n ? -log1p(n) - lbeta(n - k + 1, k + 1) + k * log(p) + (n - k) * log1p(-p) : T(-Inf)
end

binomlogpdf(n::Real, p::Real, k::Real) = binomlogpdf(promote(n, p, k)...)