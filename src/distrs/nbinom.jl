# functions related to negative binomial distribution

# R implementations
using .RFunctions:
    # nbinompdf,
    # nbinomlogpdf,
    # nbinomcdf,
    # nbinomccdf,
    # nbinomlogcdf,
    # nbinomlogccdf,
    nbinominvcdf,
    nbinominvccdf,
    nbinominvlogcdf,
    nbinominvlogccdf

function nbinompdf(r, p, x)
    exp(nbinomlogpdf(r, p, x))
end

nbinomlogpdf(n::Real, p::Real, x::Real) = nbinomlogpdf(promote(n, p, x)...)
function nbinomlogpdf(r::T, p::T, x::T) where {T <: Real}
    loggamma(r+x) - loggamma(x+1) - loggamma(r) + log1p(-p)*x + log(p)*r
end
function nbinomcdf(r, p, x)
    beta_inc(r, floor(x+1), p)[1]
end
function nbinomccdf(r, p, x)
    beta_inc(r, floor(x+1), p)[2]
end
function nbinomlogcdf(r, p, x)
    b1, b2 = beta_inc(r, floor(x+1), p)
    b1 > 1//10 ? log(b1) : log1p(-b2)
end
function nbinomlogccdf(r, p, x)
    b1, b2 = beta_inc(r, floor(x+1), p)
    b1 < 1//10 ? log1p(-b1) : log(b2)
end

# TODO: impliment https://arxiv.org/abs/2001.03953
# for inverting the incomplete beta function wrt the 2nd argument.
