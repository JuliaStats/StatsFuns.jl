# functions related to F distribution

# Julia implementations
fdistpdf(ν1::Real, ν2::Real, x::Real) = exp(fdistlogpdf(ν1, ν2, x))

fdistlogpdf(ν1::Real, ν2::Real, x::Real) = fdistlogpdf(promote(ν1, ν2, x)...)
function fdistlogpdf(ν1::T, ν2::T, x::T) where {T<:Real}
    # we ensure that `log(x)` does not error if `x < 0`
    ν1ν2 = ν1 / ν2
    y = max(x, 0)
    val = (xlogy(ν1, ν1ν2) + xlogy(ν1 - 2, y) - xlogy(ν1 + ν2, 1 + ν1ν2 * y)) / 2 - logbeta(ν1 / 2, ν2 / 2)
    return x < 0 ? oftype(val, -Inf) : val
end

for f in ("cdf", "ccdf", "logcdf", "logccdf", "invcdf", "invccdf", "invlogcdf", "invlogccdf")
    ff = Symbol("fdist"*f)
    bf = Symbol("beta"*f)
    @eval $ff(ν1::Real, ν2::Real, x::Real) = $bf(ν1/2, ν2/2, inv(1 + ν2/(ν2*max(0, x))))
end
