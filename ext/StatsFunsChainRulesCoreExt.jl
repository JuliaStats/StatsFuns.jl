module StatsFunsChainRulesCoreExt

using StatsFuns
using SpecialFunctions: digamma, erf

import ChainRulesCore

ChainRulesCore.@scalar_rule(
    betalogpdf(α::Real, β::Real, x::Number),
    @setup(z = digamma(α + β)),
    (
        log(x) + z - digamma(α),
        log1p(-x) + z - digamma(β),
        (α - 1) / x + (1 - β) / (1 - x),
    ),
)

ChainRulesCore.@scalar_rule(
    binomlogpdf(n::Real, p::Real, k::Real),
    (
        ChainRulesCore.NoTangent(),
        (k / p - n) / (1 - p),
        ChainRulesCore.NoTangent(),
    ),
)

ChainRulesCore.@scalar_rule(
    chisqlogpdf(k::Real, x::Number),
    @setup(hk = k / 2),
    (
        (log(x) - logtwo - digamma(hk)) / 2,
        (hk - 1) / x - one(hk) / 2,
    ),
)

ChainRulesCore.@scalar_rule(
    fdistlogpdf(ν1::Real, ν2::Real, x::Number),
    @setup(
        di = digamma((ν1 + ν2) / 2),
        # `r = u / (1 - u)` for the beta variate `u`, split as in `fdistlogpdf`
        ν1ν2 = ν1 / ν2,
        r = ν1ν2 * x,
        invr = inv(r),
        log1pr = isinf(r) ? log(ν1ν2) + log(x) : log1p(r),
        # `(x - 1) / (ν1 * x + ν2)`, in the form that stays finite for a large `x`
        temp1 = ν1 * x + ν2,
        a = isinf(temp1) ? (1 - inv(x)) / (ν1 + ν2 / x) : (x - 1) / temp1,
        ν2a = ν2 * a,
    ),
    (
        # `log(u)` as in `fdistlogpdf`
        ((isfinite(invr) ? -log1p(invr) : log(ν1ν2) + log(x) - log1pr) - ν2a + di - digamma(ν1 / 2)) / 2,
        (-log1pr + ν1 * a + di - digamma(ν2 / 2)) / 2,
        # `(ν1 - 2) / (2 * x) - ν1 * (ν1 + ν2) / (2 * temp1)`, combined into a single quotient
        # so that the two terms cannot cancel for a large `ν1`. `ν1 * ν2` would overflow for the
        # largest degrees of freedom, whereas `ν2 * a` is a ratio of the two
        -(ν1 * ν2a / 2 + 1) / x,
    ),
)

ChainRulesCore.@scalar_rule(
    gammalogpdf(k::Real, θ::Real, x::Number),
    @setup(
        invθ = inv(θ),
        xoθ = invθ * x,
        z = xoθ - k,
    ),
    (
        log(xoθ) - digamma(k),
        invθ * z,
        - (1 + z) / x,
    ),
)

ChainRulesCore.@scalar_rule(
    poislogpdf(λ::Number, x::Number),
    ((iszero(x) && iszero(λ) ? zero(x / λ) : x / λ) - 1, ChainRulesCore.NoTangent()),
)

ChainRulesCore.@scalar_rule(
    tdistlogpdf(ν::Real, x::Number),
    @setup(
        νp1 = ν + 1,
        xsq = x^2,
        invν = inv(ν),
        a = xsq * invν,
        b = νp1 / (ν + xsq),
    ),
    (
        (digamma(νp1 / 2) - digamma(ν / 2) + a * b - log1p(a) - invν) / 2,
        - x * b,
    ),
)

ChainRulesCore.@scalar_rule(
    owens_t(h::Real, a::Real),
    (
        normpdf(h) * erf((h * a) * invsqrt2) / -2,
        inv2π * exp(-h^2 * (1 + a^2) / 2) / (1 + a^2),
    ),
)

end # module
