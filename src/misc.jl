# Miscellaneous functions

# Logarithm of multivariate gamma function
#
# See: https://en.wikipedia.org/wiki/Multivariate_gamma_function
#
function logmvgamma(p::Int, a::Float64)
    res::Float64 = p * (p - 1.0) / 4.0 * log(pi)
    for ii in 1:p
        res += lgamma(a + (1.0 - ii) / 2.0)
    end
    return res
end

# Remainder term after Stirling's approximation to the log-gamma function
# lstirling(x) = lgamma(x) + x - (x-0.5)*log(x) - 0.5*log2π
#              = 1/(12x) - 1/(360x^3) + 1/(1260x^5) + ...
#
# Asymptotic expansion from:
#
#   Temme, N. (1996) Special functions: An introduction to the classical
#   functions of mathematical physics, Wiley, New York, ISBN: 0-471-11313-1,
#   Chapter 3.6, pp 61-65.
#
# Relative error of approximation is bounded by
#   (174611/125400 x^-19) / (1/12 x^-1 - 1/360 x^-3)
# which is < 1/2 ulp for x >= 10.0
# total numeric error appears to be < 2 ulps
#
function lstirling_asym(x::Float64)
    t = 1.0/(x*x)
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
    t = 1f0/(x*x)
    @horner(t,
             8.333333333333f-2, #  1/12 x^-1
            -2.777777777777f-3, # -1/360 x^-3
             7.936507936508f-4, #  1/1260 x^-5
            -5.952380952381f-4, # -1/1680 x^-7
             8.417508417508f-4)/x #  1/1188 x^-9
end

lstirling_asym(x::Integer) = lstirling_asym(float(x))
