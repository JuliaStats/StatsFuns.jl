module ForwardDiffExt

import StatsFuns
using ForwardDiff: ForwardDiff, Dual

# `SpecialFunctions.gamma_inc`, which every gamma tail bottoms out in, is
# defined for `Float16`, `Float32` and `Float64`, so a `Dual` needs its primal
# on the values and its tangent attached analytically.

"""
    _tail_dual(primal, derivative, k, θ, x)

Assemble the `Dual` gamma-tail result with value `primal` and `x`-derivative
`derivative`. Writing the tail as `P(k, x/θ)` gives `∂/∂θ = -(x/θ) ∂/∂x`, so
one derivative carries both tangents. `∂/∂k` has no elementary closed form,
and a nonzero shape tangent throws.
"""
function _tail_dual(
    primal::Real,
    derivative::Real,
    k::Dual{T},
    θ::Dual{T},
    x::Dual{T},
) where {T}
    if !iszero(ForwardDiff.partials(k))
        throw(ArgumentError(
            "the gamma tail functions are not differentiable with respect to " *
                "the shape parameter `k`: ∂/∂k of the regularized incomplete " *
                "gamma has no elementary closed form",
        ))
    end
    value_θ, value_x = ForwardDiff.value(θ), ForwardDiff.value(x)
    tangent = derivative * ForwardDiff.partials(x) -
        (derivative * value_x / value_θ) * ForwardDiff.partials(θ)
    return Dual{T}(primal, tangent)
end

for (f, complementary, logarithmic) in (
        (:gammacdf, false, false),
        (:gammaccdf, true, false),
        (:_gammalogcdf, false, true),
        (:_gammalogccdf, true, true),
    )
    # The log forms divide the density by their own primal; the complementary
    # forms flip its sign.
    logdensity = :(StatsFuns.gammalogpdf(value_k, value_θ, value_x))
    magnitude = if logarithmic
        :(exp($logdensity - primal))
    else
        :(exp($logdensity))
    end
    derivative = if complementary
        :(-$magnitude)
    else
        magnitude
    end
    @eval function StatsFuns.$f(k::Dual{T}, θ::Dual{T}, x::Dual{T}) where {T}
        value_k = ForwardDiff.value(k)
        value_θ = ForwardDiff.value(θ)
        value_x = ForwardDiff.value(x)
        primal = StatsFuns.$f(value_k, value_θ, value_x)
        return _tail_dual(primal, $derivative, k, θ, x)
    end
end

end # module
