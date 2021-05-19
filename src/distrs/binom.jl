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

binompdf(n::Real, p::Real, k::Real) = betapdf(k + 1, n - k + 1, p)/(n + 1)

binomlogpdf(n::Real, p::Real, k::Real) = betalogpdf(k + 1, n - k + 1, p) - log(n + 1)

binomcdf(n::Real, p::Real, k::Real) = betaccdf(k + 1, n - k, p)

binomccdf(n::Real, p::Real, k::Real) = betacdf(k + 1, n - k, p)

binomlogcdf(n::Real, p::Real, k::Real) = betalogccdf(k + 1, n - k, p)

binomlogccdf(n::Real, p::Real, k::Real) = betalogcdf(k + 1, n - k, p)
