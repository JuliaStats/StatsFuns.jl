# functions related to F distribution

# Julia implementations
fdistpdf(ν1::Real, ν2::Real, x::Real) = exp(fdistlogpdf(ν1, ν2, x))

fdistlogpdf(ν1::Real, ν2::Real, x::Real) = fdistlogpdf(promote(ν1, ν2, x)...)
function fdistlogpdf(ν1::T, ν2::T, x::T) where {T <: Real}
    a = ν1 / 2
    b = ν2 / 2
    if a isa Union{Float32, Float64} && a >= 10 && b >= 10 &&
            isfinite(a) && isfinite(b) && 0 < x < Inf
        # Older supported LogExpFunctions only specializes log1pmx for Float64.
        # Its generic Float32 subtraction loses the near-one deviance terms.
        if a isa Float32
            return Float32(_fdistlogpdf_large(Float64(a), Float64(b), Float64(x)))
        end
        return _fdistlogpdf_large(a, b, x)
    end
    lbeta = logbeta(a, b)
    # in terms of the beta variate `u = ν1 * x / (ν1 * x + ν2)` the density is
    # `u^a * (1 - u)^b / (x * beta(a, b))`; `r` below is `u / (1 - u)`, formed so that a large
    # `x` does not overflow it
    ν1ν2 = ν1 / ν2
    r = ν1ν2 * x
    invr = inv(r)
    return if x < 0
        # outside of the support, where `log(x)` would error as well
        oftype(lbeta, -Inf)
    elseif isfinite(invr)
        # `log(u) = -log1p(1 / r)` and `log(1 - u) = -log1p(r)`, symmetric in the two degrees
        # of freedom, so no two terms growing like `ν1 * log(ν1)` cancel. `r` overflows for the
        # largest `x`, where `log1p(r)` is ordinary and both summands below are positive
        blog1pr = isinf(r) ? xlogy(b, ν1ν2) + xlogy(b, x) : xlog1py(b, r)
        -xlog1py(a, invr) - blog1pr - log(x) - lbeta
    else
        # `r` is zero or subnormal: `x == 0`, `r` underflowing, or `NaN`. The textbook form
        # covers all three - `log1p(r)` keeps a tiny `r` that `log(1 + r)` would absorb, and
        # splitting `log(u)` into `log(ν1 / ν2) + log(x)` needs no precision from `r`
        xlogy(a, ν1ν2) + xlogy(a - 1, x) - xlog1py(a + b, r) - lbeta
    end
end

# Fold Stirling's formula for logbeta into the density before evaluating it.
# With p = a / (a + b), q = b / (a + b), d = p*x + q, t = q*(x - 1)/d
# and z = -p*(x - 1)/d, a*t + b*z = 0. Thus the large cancelling terms
# reduce to a*log1pmx(t) + b*log1pmx(z), both nonpositive, even near x = 1.
function _fdistlogpdf_large(a, b, x)
    s = a + b # finite: a and b are halves of finite degrees of freedom
    p, q = a / s, b / s
    d = p * x + q
    t = q * (x - 1) / d
    z = -p * (x - 1) / d
    # Near -1 the residual loses its positive complement. Evaluate the log
    # ratio directly there, rather than forming 1 + t or 1 + z.
    da = t > -1 // 2 ? log1pmx(t) : log(x) - log(d) - t
    db = z > -1 // 2 ? log1pmx(z) : -log(d) - z
    small, large = minmax(a, b)
    logharmonic = log(small) - log1p(small / large) # log(a*b/(a+b))
    return a * da + b * db - log(x) + (logharmonic - oftype(a, log2π)) / 2 -
        lstirling_asym(a) - lstirling_asym(b) + lstirling_asym(s)
end

for f in ("cdf", "ccdf", "logcdf", "logccdf")
    ff = Symbol("fdist" * f)
    bf = Symbol("beta" * f)
    @eval function $ff(ν1::T, ν2::T, x::T) where {T <: Real}
        # the beta variate `u = r / (1 + r)`, clamped to the support. Here an overflowing
        # `ν1 * x` only saturates `u` to 1, so the ratio is formed that way round
        r = (ν1 * max(0, x)) / ν2
        u = r > 1 ? inv(1 + inv(r)) : r / (1 + r)
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
