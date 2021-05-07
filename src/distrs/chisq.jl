# functions related to chi-square distribution

# R implementations
using .RFunctions:
    # chisqpdf,
    # chisqlogpdf,
    # chisqcdf,
    # chisqccdf,
    chisqlogcdf,
    # chisqlogccdf,
    # chisqinvcdf,
    # chisqinvccdf,
    chisqinvlogcdf,
    chisqinvlogccdf

# Julia implementations
# promotion ensures that we do forward e.g. `chisqpdf(::Int, ::Float32)` to
# `gammapdf(::Float32, ::Int, ::Float32)` but not `gammapdf(::Float64, ::Int, ::Float32)`
chisqpdf(k::Real, x::Real) = chisqpdf(promote(k, x)...)
chisqpdf(k::T, x::T) where {T<:Real} = gammapdf(k / 2, 2, x)

chisqlogpdf(k::Real, x::Real) = chisqlogpdf(promote(k, x)...)
chisqlogpdf(k::T, x::T) where {T<:Real} = gammalogpdf(k / 2, 2, x)

chisqcdf(k::Float64, x::Float64) = first(gamma_inc(k/2, x/2, 0))
chisqcdf(k::Real, x::Real) = chisqcdf(promote(float(k), x)...)
chisqcdf(k::T, x::T) where T = throw(MethodError(chisqcdf, (k, x)))

chisqccdf(k::Float64, x::Float64) = last(gamma_inc(k/2, x/2, 0))
chisqccdf(k::Real, x::Real) = chisqccdf(promote(float(k), x)...)
chisqccdf(k::T, x::T) where T = throw(MethodError(chisqccdf, (k, x)))

chisqlogccdf(k::Float64, x::Float64) = loggamma(k/2, x/2) - loggamma(k/2)
chisqlogccdf(k::Real, x::Real) = chisqlogccdf(promote(float(k), x)...)
chisqlogccdf(k::T, x::T) where T = throw(MethodError(chisqlogccdf, (k, x)))

chisqinvcdf(k::Float64, p::Float64) = 2*gamma_inc_inv(k/2, p, 1 - p)
chisqinvcdf(k::Real, p::Real) = chisqinvcdf(promote(float(k), p)...)
chisqinvcdf(k::T, p::T) where T = throw(MethodError(chisqinvcdf, (k, p)))

chisqinvccdf(k::Float64, p::Float64) = 2*gamma_inc_inv(k/2, 1 - p, p)
chisqinvccdf(k::Real, p::Real) = chisqinvccdf(promote(float(k), p)...)
chisqinvccdf(k::T, p::T) where T = throw(MethodError(chisqinvccdf, (k, p)))
