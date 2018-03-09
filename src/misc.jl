# Miscellaneous functions

"""
    logmvgamma(p::Int, a::Real)

Return the logarithm of [multivariate gamma function](https://en.wikipedia.org/wiki/Multivariate_gamma_function)
"""
function logmvgamma(p::Int, a::Real)
    # NOTE: one(a) factors are here to prevent unnecessary promotion of Float32
    res = p * (p - 1) * log(pi * one(a)) / 4
    for ii in 1:p
        res += lgamma(a + (1 - ii) * one(a)/ 2)
    end
    return res
end

"""
    lstirling_asym(x)

Return the remainder term after Stirling's approximation to `lgamma(x)`.

[Stirling's approximation](https://en.wikipedia.org/wiki/Stirling%27s_approximation)
to the log-gamma (`lgamma`) function also applies to `log(factorial(k))`  through
the relationship `log(factorial(k)) = lgamma(k + 1)`.  Stirling's approximation is

    lgamma(x) ≈ x*log(x) - x + log(2π/x)/2 = log(x)*(x-0.5) + log2π/2 -x

Thus

    lstirling_asym(x) = lgamma(x) + x - (x-0.5)*log(x) - 0.5*log2π
                      = 1/(12x) - 1/(360x^3) + 1/(1260x^5) + ...

Why bother with this when `lgamma` can be evaluated directly for most numeric types?
In some formulas it is advantageous to evaluate the Stirling approximation and
the remainder to take advantage of explicit cancelation of some terms in
the Stirling approximation.

Asymptotic expansion from:

   Temme, N. (1996) Special functions: An introduction to the classical
   functions of mathematical physics, Wiley, New York, ISBN: 0-471-11313-1,
   Chapter 3.6, pp 61-65.

Relative error of approximation is bounded by
    (174611/125400 x^-19) / (1/12 x^-1 - 1/360 x^-3)
which is < 1/2 ulp for x >= 10.0
total numeric error appears to be < 2 ulps
"""
function lstirling_asym end

lstirling_asym(x::BigFloat) = lgamma(x) + x - log(x)*(x - big(0.5)) - log2π/big(2)

lstirling_asym(x::Integer) = lstirling_asym(float(x))

const lstirlingF64 = Float64[lstirling_asym(k) for k in big(1):big(64)]
const lstirlingF32 = Float64[lstirling_asym(k) for k in big(1):big(40)]

function lstirling_asym(x::Float64)
    isinteger(x) && (0 < x ≤ length(lstirlingF64)) && return lstirlingF64[Int(x)]
    t = inv(abs2(x))
    @horner(t,
             8.33333333333333333e-2, #  1/12 x^-1
            -2.77777777777777778e-3, # -1/360 x^-3
             7.93650793650793651e-4, #  1/1260 x^-5
            -5.95238095238095238e-4, # -1/1680 x^-7
             8.41750841750841751e-4, #  1/1188 x^-9
            -1.91752691752691753e-3, # -691/360360 x^-11
             6.41025641025641026e-3, #  1/156 x^-13
            -2.95506535947712418e-2, # -3617/122400 x^-15
             1.79644372368830573e-1)/x #  43867/244188 x^-17
end

function lstirling_asym(x::Float32)
    isinteger(x) && (0 < x ≤ length(lstirlingF32)) && return lstirlingF32[Int(x)]
    t = inv(abs2(x))
    @horner(t,
             8.333333333333f-2, #  1/12 x^-1
            -2.777777777777f-3, # -1/360 x^-3
             7.936507936508f-4, #  1/1260 x^-5
            -5.952380952381f-4, # -1/1680 x^-7
             8.417508417508f-4)/x #  1/1188 x^-9
end
