# functions related to student's T distribution

# R implementations
using .RFunctions:
    tdistinvlogcdf,
    tdistinvlogccdf

# Julia implementations
tdistpdf(ν::Real, x::Real) = exp(tdistlogpdf(ν, x))

tdistlogpdf(ν::Real, x::Real) = tdistlogpdf(promote(ν, x)...)
function tdistlogpdf(ν::T, x::T) where {T<:Real}
    isinf(ν) && return normlogpdf(x)
    νp12 = (ν + 1) / 2
    return loggamma(νp12) - (logπ + log(ν)) / 2 - loggamma(ν / 2) - νp12 * log1p(x^2 / ν)
end
function tdistcdf(ν::T, x::T) where {T<:Real}
    if x < 0 && ν < (0x1p-26 * x)^2
       return -log(abs(x))*ν + log(ν)*muladd(ν, T(.5), -1.) - logbeta(ν/2, T(.5))
    end
    q = T(0.5)*beta_inc(ν/2, T(.5), ν/muladd(x, x, ν))[1]
    return ifelse(q > 0, 1 - q, q)
end
function tdistlogcdf(ν::T, x::T) where {T<:Real}
    if x < 0 && ν < (0x1p-26 * x)^2
       return -log(abs(x))*ν + log(ν)*muladd(ν, T(.5), -1) - logbeta(ν/2, T(.5))
    end
    q = T(.5)*beta_inc(ν/2, T(.5), ν/muladd(x, x, r))[1]
    return q > 0 ? log1p(-q) : log(q)
end
tdistccdf(ν, x) = tdistcdf(ν, -x)
tdistlogccdf(ν, x) = tdistlogcdf(ν, -x)

function tdistinvcdf(ν, x::Real)
    if x > .5
        return -tdistinvcdf(d, 1-x)
    end
    return -sqrt(ν * (inv(beta_inc_inv(ν/2, 1/2, 2*x)[1]) - 1))
end
tdistinvccdf(ν, x::Real) = -tdistinvcdf(ν, x)
