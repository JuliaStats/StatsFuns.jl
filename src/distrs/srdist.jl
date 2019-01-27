# functions related to studentized range distribution

using ForwardDiff: derivative
using QuadGK: quadgk

import .RFunctions:
    srdistcdf,
    srdistccdf,
    srdistlogcdf,
    srdistlogccdf,
    srdistinvcdf,
    srdistinvccdf,
    srdistinvlogcdf,
    srdistinvlogccdf

# pdf for numbers with generic types
function srdistpdf(ν::Real, k::Real, x::Number)
    Φ(x) = (1+erf(x / √2)) / 2
    ϕ(x) = derivative(Φ, x)
    function outer(y)
        function inner(u)
            return ϕ(u) * ϕ(u - x*y) * (Φ(u) - Φ(u - x*y))^(k-2)
        end
        inner_part = quadgk(inner, -Inf, Inf)[1]
        return inner_part * y^ν * ϕ(y*√ν)
    end
    integral = quadgk(outer, 0.0, Inf)[1]
    return integral * (√(2π) * k * (k-1) * ν^(ν/2)) / (gamma(ν/2) * 2^(ν/2 - 1))
end

# logpdf for numbers with generic types
srdistlogpdf(ν::Real, k::Real, x::Number) = log(srdistpdf(ν, k, x))
