# functions related to chi-square distribution

# Just use the Gamma definitions
for f in ("pdf", "logpdf", "cdf", "ccdf", "logcdf", "logccdf", "invcdf", "invccdf", "invlogcdf", "invlogccdf")
    _chisqf = Symbol("chisq"*f)
    _gammaf = Symbol("gamma"*f)
    @eval $(_chisqf)(k::Real, x::Real) = $(_gammaf)(k/2, 2, x)
end
