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
binomlogpdf(n::Real, p::Real, k::Real) = (isinteger(k) & (zero(k) <= k <= n)) ? convert(typeof(float(p)), loggamma(n + 1) - loggamma(k + 1) - loggamma(n - k + 1) + xlogy(k, p) + xlogy(n - k, 1 - p)) : convert(typeof(float(p)), -Inf)
