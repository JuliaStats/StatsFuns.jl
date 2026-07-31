# functions related to F distribution

# Julia implementations
fdistpdf(ν1::Real, ν2::Real, x::Real) = exp(fdistlogpdf(ν1, ν2, x))

fdistlogpdf(ν1::Real, ν2::Real, x::Real) = fdistlogpdf(promote(ν1, ν2, x)...)
function fdistlogpdf(ν1::T, ν2::T, x::T) where {T <: Real}
    lbeta = logbeta(ν1 / 2, ν2 / 2)
    return if x > 0
        # in terms of the beta variate `u = ν1 * x / (ν1 * x + ν2)`, which is symmetric in
        # `ν1` and `ν2`: neither is absorbed by an explicit `1 +`, nor cancels the other
        -xlog1py(ν1 / 2, ν2 / (ν1 * x)) - xlog1py(ν2 / 2, ν1 * x / ν2) - log(x) - lbeta
    elseif x < 0
        oftype(lbeta, -Inf)
    else
        # at zero the density behaves like `x^(ν1 / 2 - 1)`; NaN propagates through `xlogy`
        (xlogy(ν1, ν1 / ν2) + xlogy(ν1 - 2, x)) / 2 - lbeta
    end
end

for f in ("cdf", "ccdf", "logcdf", "logccdf")
    ff = Symbol("fdist" * f)
    bf = Symbol("beta" * f)
    @eval $ff(ν1::T, ν2::T, x::T) where {T <: Real} = $bf(ν1 / 2, ν2 / 2, inv(1 + ν2 / (ν1 * max(0, x))))
    @eval $ff(ν1::Real, ν2::Real, x::Real) = $ff(promote(ν1, ν2, x)...)
end
for f in ("invcdf", "invccdf", "invlogcdf", "invlogccdf")
    ff = Symbol("fdist" * f)
    bf = Symbol("beta" * f)
    @eval function $ff(ν1::T, ν2::T, y::T) where {T <: Real}
        x = $bf(ν1 / 2, ν2 / 2, y)
        return x / (1 - x) * ν2 / ν1
    end
    @eval $ff(ν1::Real, ν2::Real, y::Real) = $ff(promote(ν1, ν2, y)...)
end
