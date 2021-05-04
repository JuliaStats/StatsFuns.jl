# functions related to chi-square distribution

import .RFunctions:
    # chisqpdf,
    # chisqlogpdf,
    chisqcdf,
    chisqccdf,
    chisqlogcdf,
    chisqlogccdf,
    chisqinvcdf,
    chisqinvccdf,
    chisqinvlogcdf,
    chisqinvlogccdf

# pdf for numbers with generic types
function chisqpdf(k::T, x::T) where T<:Real
  hk = k / 2  # half k
  return 1 / (2^(hk) * gamma(hk)) * x^(hk - 1) * exp(-x / 2)
end
chisqpdf(k::Real, x::Real) = chisqpdf(promote(k, x)...)

# logpdf for numbers with generic types
function chisqlogpdf(k::T, x::T) where T<:Real
  hk = k / 2  # half k
  return -hk * log(oftype(hk, 2)) - loggamma(hk) + (hk - 1) * log(x) - x / 2
end
chisqlogpdf(k::Real, x::Real) = chisqlogpdf(promote(k, x)...)
