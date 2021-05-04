# functions related to chi-square distribution

# R implementations
using .RFunctions:
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

# Julia implementations
# promotion ensures that we do forward e.g. `chisqpdf(::Int, ::Float32)` to
# `gammapdf(::Float32, ::Int, ::Float32)` but not `gammapdf(::Float64, ::Int, ::Float32)`
chisqpdf(k::Real, x::Real) = chisqpdf(promote(k, x)...)
chisqpdf(k::T, x::T) where {T<:Real} = gammapdf(k / 2, 2, x)

chisqlogpdf(k::Real, x::Real) = chisqlogpdf(promote(k, x)...)
chisqlogpdf(k::T, x::T) where {T<:Real} = gammalogpdf(k / 2, 2, x)
