# functions related to binomial distribution

# R implementations
using .RFunctions:
    # binompdf,
    # binomlogpdf,
    # binomcdf,
    # binomccdf,
    # binomlogcdf,
    # binomlogccdf,
    binominvcdf,
    binominvccdf,
    binominvlogcdf,
    binominvlogccdf

# Julia implementations
binompdf(n::Real, p::Real, k::Real) = exp(binomlogpdf(n, p, k))

binomlogpdf(n::Real, p::Real, k::Real) = binomlogpdf(promote(n, p, k)...)
function binomlogpdf(n::T, p::T, k::T) where {T<:Real}
    m = clamp(k, 0, n)
    val = min(0, betalogpdf(m + 1, n - m + 1, p) - log(n + 1))
    return 0 <= k <= n && isinteger(k) ? val : oftype(val, -Inf)
end

binomcdf(n::Real, p::Real, k::Real) = betaccdf(k + 1, n - k, p)

binomccdf(n::Real, p::Real, k::Real) = betacdf(k + 1, n - k, p)

binomlogcdf(n::Real, p::Real, k::Real) = betalogccdf(k + 1, n - k, p)

binomlogccdf(n::Real, p::Real, k::Real) = betalogcdf(k + 1, n - k, p)
