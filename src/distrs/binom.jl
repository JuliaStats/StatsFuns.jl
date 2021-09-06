# functions related to binomial distribution

# R implementations
# For pdf and logpdf we use the Julia implementation
using .RFunctions:
    binomcdf,
    binomccdf,
    binomlogcdf,
    binomlogccdf,
    binominvcdf,
    binominvccdf,
    binominvlogcdf,
    binominvlogccdf


# Julia implementations
binompdf(n::Real, p::Real, k::Real) = exp(binomlogpdf(n, p, k))

binomlogpdf(n::Real, p::Real, k::Real) = binomlogpdf(promote(n, p, k)...)
function binomlogpdf(n::T, p::T, k::T) where {T<:Real}
    val = betalogpdf(k + 1, n - k + 1, p) - log(n + 1)
    return 0 <= k <= n && isinteger(n) && isinteger(k) ? val : oftype(val, -Inf)
end
