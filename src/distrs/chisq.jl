# functions related to chi-square distribution

# Just use the Gamma definitions
for f in ("pdf", "logpdf", "cdf", "ccdf", "logcdf", "logccdf", "invcdf", "invccdf", "invlogcdf", "invlogccdf")
    _chisqf = Symbol("chisq"*f)
    _gammaf = Symbol("gamma"*f)
    @eval begin
        $(_chisqf)(k::Real, x::Real) = $(_chisqf)(promote(k, x)...)
        $(_chisqf)(k::T, x::T) where {T<:Real} = $(_gammaf)(k/2, 2, x)
    end
end
