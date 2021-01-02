# functions related to chi-square distribution

import .RFunctions:
    chisqpdf,
    chisqlogpdf,
    chisqcdf,
    chisqccdf,
    chisqlogcdf,
    chisqlogccdf,
    chisqinvcdf,
    chisqinvccdf,
    chisqinvlogcdf,
    chisqinvlogccdf

# pdf for numbers with generic types
function chisqpdf(k::Real, x::Number)
  hk = k / 2  # half k
  1 / (2^(hk) * gamma(hk)) * x^(hk - 1) * exp(-x / 2)
end

# logpdf for numbers with generic types
function chisqlogpdf(k::Real, x::Number)
  hk = k / 2  # half k
  -hk * logtwo - loggamma(hk) + (hk - 1) * log(x) - x / 2
end

# ChainRules adjoint
ChainRulesCore.@scalar_rule(
    chisqlogpdf(k::Real, x::Number),
    @setup(hk = k / 2),
    (
        (log(x) - logtwo - digamma(hk)) / 2,
        (hk - 1) / x - one(hk) / 2,
    ),
)
