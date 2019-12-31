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
function binomlogpdf(n::Real, p::Real, k::Real)
    if isinteger(k) & (zero(k) <= k <= n)
        x = loggamma(n + 1) - loggamma(k + 1) - loggamma(n - k + 1) +
            xlogy(k, p) + xlogy(n - k, 1 - p)
    else
        x = -Inf
    end

    return convert(typeof(float(p)), x)
end
