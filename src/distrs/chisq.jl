# functions related to chi-square distribution

# Just use the Gamma definitions
for f in ("pdf", "logpdf", "logupdf", "cdf", "ccdf", "logcdf", "logccdf", "invcdf", "invccdf", "invlogcdf", "invlogccdf")
    _chisqf = Symbol("chisq" * f)
    _gammaf = Symbol("gamma" * f)
    @eval begin
        $(_chisqf)(k::Real, x::Real) = $(_chisqf)(promote(k, x)...)
        $(_chisqf)(k::T, x::T) where {T <: Real} = $(_gammaf)(k / 2, 2, x)
    end
end

chisqlogulikelihood(k::Real, x::Real) = chisqlogulikelihood(promote(k, x)...)
function chisqlogulikelihood(k::T, x::T) where {T <: Real}
    y = max(x, 0)
    k2 = k / 2
    val = xlogy(k2, x / 2) - loggamma(k2)
    return x < 0 ? oftype(val, -Inf) : val
end
