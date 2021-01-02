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
binomlogpdf(n::Real, p::Real, k::Real) = -log1p(n) - logbeta(n - k + 1, k + 1) + k * log(p) + (n - k) * log1p(-p)

# ChainRules adjoint
ChainRulesCore.@scalar_rule(
    binomlogpdf(n::Real, p::Real, k::Real),
    @setup(z = digamma(n - k + 1)),
    (
        digamma(n + 2) - z + log1p(-p) - 1 / (1 + n),
        (k / p - n) / (1 - p),
        z - digamma(k + 1) + logit(p),
    ),
)

