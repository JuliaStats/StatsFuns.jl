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
chisqpdf(k::Real, x::Real) = exp(chisqlogpdf(k, x))

# logpdf for numbers with generic types
function chisqlogpdf(k::T, x::T) where {T <: Real}
  isinteger(k) && k > 0 || throw(ArgumentError("k must be a positive integer, got $k"))
  x ≥ 0 || return T(-Inf)
  hk = k / 2  # half k
  -hk * logtwo - lgamma(hk) + (hk - 1) * log(x) - x / 2
end

chisqlogpdf(k::Real, x::Real) = chisqlogpdf(promote(k, float(x))...)
