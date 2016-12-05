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
  hk = 0.5 * k  # hlaf k
  1 / (2^(hk) * gamma(hk)) * x^(hk - 1) * exp(-0.5 * x)
end

# logpdf for numbers with generic types
function chisqlogpdf(k::Real, x::Number)
  hk = 0.5 * k  # hlaf k
  -hk * log(2) - lgamma(hk) + (hk - 1) * log(x) - 0.5 * x
end
