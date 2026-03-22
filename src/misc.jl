# Miscellaneous functions

"""
    logmvgamma(p::Int, a::Real)

Return the logarithm of [multivariate gamma function](https://en.wikipedia.org/wiki/Multivariate_gamma_function) ([DLMF 35.3.1](https://dlmf.nist.gov/35.3.1)).
"""
function logmvgamma(p::Int, a::Real)
    # NOTE: one(a) factors are here to prevent unnecessary promotion of Float32
    res = p * (p - 1) * (logπ * one(a)) / 4
    for ii in 1:p
        res += loggamma(a + (1 - ii) * one(a) / 2)
    end
    return res
end

"""
    logmvbeta(p::Int, a::Real, b::Real)

Return the logarithm of the multivariate beta function ([DLMF 35.3.7](https://dlmf.nist.gov/35.3#E7)).
"""
logmvbeta(p::Int, a::T, b::T) where {T <: Real} = logmvgamma(p, a) + logmvgamma(p, b) - logmvgamma(p, a + b)
logmvbeta(p::Int, a::Real, b::Real) = logmvbeta(p, promote(a, b)...)

"""
    lstirling_asym(x)

The remainder term after
[Stirling's approximation](https://en.wikipedia.org/wiki/Stirling%27s_approximation)
to [`loggamma`](@ref):

```math
\\log \\Gamma(x) \\approx x \\log(x) - x + \\log(2\\pi/x)/2 = \\log(x)*(x-1/2) + \\log(2\\pi)/2 - x
```

In Julia syntax, this means:

    lstirling_asym(x) = loggamma(x) + x - (x-0.5)*log(x) - 0.5*log(2π)

For sufficiently large `x`, this can be approximated using the asymptotic
_Stirling's series_ ([DLMF 5.11.1](https://dlmf.nist.gov/5.11.1)):

```math
\\frac{1}{12x} - \\frac{1}{360x^3} + \\frac{1}{1260x^5} - \\frac{1}{1680x^7} + \\ldots
```

The truncation error is bounded by the first omitted term, and is of the same sign.

Relative error of approximation is bounded by
    (174611/125400 x^-19) / (1/12 x^-1 - 1/360 x^-3)
which is < 1/2 ulp for x >= 10.0, and total numeric error appears to be < 2 ulps

# References

* Temme, N. (1996) Special functions: An introduction to the classical functions of
   mathematical physics, Wiley, New York, ISBN: 0-471-11313-1, Chapter 3.6, pp 61-65.
* Weisstein, Eric W. ["Stirling's Series."](http://mathworld.wolfram.com/StirlingsSeries.html).
  MathWorld.
* [OEIS A046968](http://oeis.org/A046968) and [OEIS A046969](http://oeis.org/A046969)
  for the series coefficients
"""
function lstirling_asym end

lstirling_asym(x::BigFloat) = loggamma(x) + x - log(x) * (x - big(0.5)) - log2π / big(2)

lstirling_asym(x::Integer) = lstirling_asym(float(x))

const lstirlingF64 = Float64[lstirling_asym(k) for k in big(1):big(64)]
const lstirlingF32 = [Float32(lstirlingF64[i]) for i in 1:40]

function lstirling_asym(x::Float64)
    isinteger(x) && (0 < x ≤ length(lstirlingF64)) && return lstirlingF64[Int(x)]
    t = inv(abs2(x))
    return @horner(
        t,
        8.33333333333333333e-2, #  1/12 x^-1
        -2.77777777777777778e-3, # -1/360 x^-3
        7.93650793650793651e-4, #  1/1260 x^-5
        -5.95238095238095238e-4, # -1/1680 x^-7
        8.41750841750841751e-4, #  1/1188 x^-9
        -1.91752691752691753e-3, # -691/360360 x^-11
        6.41025641025641026e-3, #  1/156 x^-13
        -2.95506535947712418e-2, # -3617/122400 x^-15
        1.79644372368830573e-1
    ) / x #  43867/244188 x^-17
end

function lstirling_asym(x::Float32)
    isinteger(x) && (0 < x ≤ length(lstirlingF32)) && return lstirlingF32[Int(x)]
    t = inv(abs2(x))
    return @horner(
        t,
        8.333333333333f-2, #  1/12 x^-1
        -2.777777777777f-3, # -1/360 x^-3
        7.936507936508f-4, #  1/1260 x^-5
        -5.952380952381f-4, # -1/1680 x^-7
        8.417508417508f-4
    ) / x #  1/1188 x^-9
end

"""
    logfbit(x)

Stirling error term for log-factorial:

    logfbit(x) = log(x!) - log(√2π) + (x+1) - (x+0.5)*log(x+1)

Equivalent to `lstirling_asym(x + 1)`.
"""
logfbit(x) = lstirling_asym(x + one(x))

"""
    lfbaccdif1(a, b)

Accurate computation of `logfbit(b) - logfbit(a + b)`.
Uses a polynomial expansion for `b ≥ 8` that avoids cancellation.
Based on VBA code by Ian Smith.
"""
function lfbaccdif1(a::Float64, b::Float64)
    if a < 0
        return -lfbaccdif1(-a, b + a)
    end
    if b >= 8
        y1 = b + 1.0
        y2 = inv(y1 * y1)
        x1 = a + b + 1.0
        x2 = inv(x1 * x1)

        # Initialize with innermost tuned coefficient (lfbc9)
        x3 = x2 * 1.6769380337122674863
        y3 = y2 * 1.6769380337122674863
        acc = x2 * (a * (x1 + y1) * y3)

        # Unroll from lfbc8 down to lfbc2
        for c in (0.35068485511628418514, 1 / 13, 691 / 30030, 1 / 99, 1 / 140, 1 / 105, 1 / 30)
            x3 = x2 * (c - x3)
            y3 = y2 * (c - y3)
            acc = x2 * (a * (x1 + y1) * y3 - acc)
        end

        return (a * (1.0 - y3) - y1 * acc) / (12.0 * x1 * y1)
    else
        return logfbit(b) - logfbit(a + b)
    end
end
lfbaccdif1(a::Real, b::Real) = lfbaccdif1(Float64(a), Float64(b))

"""
    ab_minus_cd(a, b, c, d)

Accurate computation of `a * b - c * d` using FMA.
"""
function ab_minus_cd(a::Float64, b::Float64, c::Float64, d::Float64)
    w = c * d
    return fma(a, b, -w) - fma(c, d, -w)
end
ab_minus_cd(a::Real, b::Real, c::Real, d::Real) = ab_minus_cd(Float64(a), Float64(b), Float64(c), Float64(d))

# Noncentral distribution helpers (based on VBA code by Ian Smith)

const _minLog1Value = -0.79149064

"""
    _logdif(pr, prob)

Accurate computation of `log(pr / prob)`. Uses `log1pmx`-based computation
when `pr` is close to `prob` to avoid cancellation.
Based on VBA `logdif` by Ian Smith.
"""
function _logdif(pr::Float64, prob::Float64)
    temp = (pr - prob) / prob
    if abs(temp) >= 0.5
        return log(pr / prob)
    else
        return log1p(temp) # log(1 + temp) = log(pr/prob)
    end
end

"""
    _poisson_term(i, n, diffFromMean, logAdd)

High-precision Poisson PMF: probability that a Poisson variate with mean `n`
has value `i`, where `diffFromMean = n - i`. The result is multiplied by
`exp(logAdd)`. Uses Stirling corrections for accuracy.
Based on VBA `poissonTerm` by Ian Smith.
"""
function _poisson_term(i::Float64, n::Float64, diffFromMean::Float64, logAdd::Float64)
    if (i <= -1.0) || (n < 0.0)
        if i == 0.0
            return exp(logAdd)
        else
            return 0.0
        end
    elseif (i < 0.0) && (n == 0.0)
        return NaN
    else
        c3 = i
        c2 = c3 + 1.0
        c1 = (diffFromMean - 1.0) / c2
        if c1 < _minLog1Value
            if i == 0.0
                logpoissonTerm = -n
                return exp(logpoissonTerm + logAdd)
            else
                logpoissonTerm = (c3 * log(n / c2) - (diffFromMean - 1.0)) - logfbit(c3)
                r = exp(logpoissonTerm + logAdd)
                (isfinite(r)) || return 0.0
                return r / sqrt(c2) * Float64(invsqrt2π)
            end
        else
            logpoissonTerm = c3 * log1pmx(c1) - c1 - logfbit(c3)
            return exp(logpoissonTerm + logAdd) / sqrt(c2) * Float64(invsqrt2π)
        end
    end
end

"""
    _binomial_term(i, j, p, q, diffFromMean, logAdd)

High-precision binomial PMF: probability that a binomial variate with sample
size `i + j` and event probability `p` (where `q = 1 - p`) has value `i`,
where `diffFromMean = (i + j) * p - i`. The result is multiplied by `exp(logAdd)`.
Uses Stirling corrections for accuracy.
Based on VBA `binomialTerm` by Ian Smith.
"""
function _binomial_term(i::Float64, j::Float64, p::Float64, q::Float64, diffFromMean::Float64, logAdd::Float64)
    if (i == 0.0) && (j <= 0.0)
        return exp(logAdd)
    elseif (i <= -1.0) || (j < 0.0)
        return 0.0
    else
        if p < q
            c2 = i
            c3 = j
            ps = p
            dfm = diffFromMean
        else
            c3 = i
            c2 = j
            ps = q
            dfm = -diffFromMean
        end
        c5 = (dfm - (1.0 - ps)) / (c2 + 1.0)
        c6 = -(dfm + ps) / (c3 + 1.0)
        if c5 < _minLog1Value
            if c2 == 0.0
                logbinomialTerm = c3 * log1p(-ps)
                return exp(logbinomialTerm + logAdd)
            elseif (ps == 0.0) && (c2 > 0.0)
                return 0.0
            else
                c1 = (i + 1.0) + j
                c4 = lfbaccdif1(j, i) + logfbit(j)
                logbinomialTerm = c2 * (log((ps * c1) / (c2 + 1.0)) - c5) - c5 + c3 * log1pmx(c6) - c6 - c4
                return exp(logbinomialTerm + logAdd) * sqrt(c1 / ((c2 + 1.0) * (c3 + 1.0))) * Float64(invsqrt2π)
            end
        else
            c4 = lfbaccdif1(j, i) + logfbit(j)
            logbinomialTerm = (c2 * log1pmx(c5) - c5) + (c3 * log1pmx(c6) - c6) - c4
            return exp(logbinomialTerm + logAdd) * sqrt((1.0 + j / (i + 1.0)) / (j + 1.0)) * Float64(invsqrt2π)
        end
    end
end
