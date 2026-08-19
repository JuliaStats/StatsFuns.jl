# functions related to F distribution

# Julia implementations
fdistpdf(ν1::Real, ν2::Real, x::Real) = exp(fdistlogpdf(ν1, ν2, x))

fdistlogpdf(ν1::Real, ν2::Real, x::Real) = fdistlogpdf(promote(ν1, ν2, x)...)
function fdistlogpdf(ν1::T, ν2::T, x::T) where {T <: Real}
    a = ν1 / 2
    b = ν2 / 2
    lbeta = logbeta(a, b)
    # in terms of the beta variate `u = ν1 * x / (ν1 * x + ν2)` the density is
    # `u^a * (1 - u)^b / (x * beta(a, b))`; `r` below is `u / (1 - u)`
    r = ν1 * x / ν2
    return if x < 0
        # outside of the support, where `log(x)` would error as well
        oftype(lbeta, -Inf)
    elseif r > 1
        # `u > 1/2`: `log(u) = -log1p(1 / r)` and `log(1 - u) = -log1p(r)`, which is symmetric
        # in the two degrees of freedom, so no two terms growing like `ν1 * log(ν1)` cancel
        -xlog1py(a, inv(r)) - xlog1py(b, r) - log(x) - lbeta
    else
        # `u <= 1/2`: the textbook form, but with `log1p(r)` instead of `log(1 + r)`, which
        # would absorb a tiny `r`. Splitting `log(u)` into `log(ν1 / ν2) + log(x)` also keeps
        # it accurate if `r` underflows. Covers `x == 0`, where the density behaves like
        # `x^(a - 1)`, and `NaN`, which propagates through `xlogy`.
        xlogy(a, ν1 / ν2) + xlogy(a - 1, x) - xlog1py(a + b, r) - lbeta
    end
end

for f in ("cdf", "ccdf", "logcdf", "logccdf")
    ff = Symbol("fdist" * f)
    bf = Symbol("beta" * f)
    @eval function $ff(ν1::T, ν2::T, x::T) where {T <: Real}
        # the beta variate `u = y / (y + ν2)`, clamped to the support. `ν2 / y` overflows for
        # small `y` and `y / (y + ν2)` is `NaN` for `y = Inf`, so each form is only used where
        # it holds up
        y = ν1 * max(0, x)
        u = y > ν2 ? inv(1 + ν2 / y) : y / (y + ν2)
        return $bf(ν1 / 2, ν2 / 2, u)
    end
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
