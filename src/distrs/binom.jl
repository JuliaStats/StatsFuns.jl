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
binomlogpdf(n::Real, p::Real, k::Real) = begin
    -log(n + 1) - lbeta(n - k + 1, k + 1) + k * log(p) + (n - k) * log(1 - p)
end
