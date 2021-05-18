# functions related to binomial distribution

import .RFunctions:
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

# pdf for numbers with generic types
binompdf(n::Real, p::Real, k::Real) = exp(binomlogpdf(n, p, k))

# logpdf for numbers with generic types
binomlogpdf(n::T, p::T, k::T) where T<:Real = -log1p(n) - logbeta(n - k + 1, k + 1) + k * log(p) + (n - k) * log1p(-p)
binomlogpdf(n::Real, p::Real, k::Real) = binomlogpdf(promote(n, p, k)...)

binomcdf(n::Real, p::Real, k::Real) = betaccdf(k + 1, n - k, p)

binomccdf(n::Real, p::Real, k::Real) = betacdf(k + 1, n - k, p)

binomlogcdf(n::Real, p::Real, k::Real) = betalogccdf(k + 1, n - k, p)

binomlogccdf(n::Real, p::Real, k::Real) = betalogcdf(k + 1, n - k, p)
