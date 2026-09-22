module ForwardDiffExt

import StatsFuns
using ForwardDiff: ForwardDiff

# `SpecialFunctions.gamma_inc`, which every gamma tail bottoms out in, is defined
# for `Float16`, `Float32` and `Float64` only, so a `Dual` argument needs its
# primal on the values and its tangent attached analytically.

"""
    _tail_partials(primal, density, θ, x, complementary, logarithmic)

The `θ` and `x` partials of a gamma tail whose value is `primal`.

Writing the tail as `P(k, x/θ)` gives `∂/∂θ = -(x/θ) ∂/∂x`, and `∂/∂x` is the
density, negated for the complementary tail and divided by the value for the log
forms.
"""
@inline function _tail_partials(
    primal::Real,
    logdensity::Real,
    θ::Real,
    x::Real,
    complementary::Bool,
    logarithmic::Bool,
)
    magnitude = logarithmic ? exp(logdensity - primal) : exp(logdensity)
    derivative_x = complementary ? -magnitude : magnitude
    return -(x / θ) * derivative_x, derivative_x
end

# `∂/∂k` of the regularized incomplete gamma has no elementary closed form, so a
# dual shape parameter is an error rather than a silently dropped tangent.
@noinline _shape_is_not_differentiable() = throw(ArgumentError(
    "the gamma tail functions are not differentiable with respect to the " *
        "shape parameter",
))

for (f, complementary, logarithmic) in (
        (:gammacdf, false, false),
        (:gammaccdf, true, false),
        (:gammalogcdf, false, true),
        (:gammalogccdf, true, true),
    )
    @eval begin
        @inline function _dual_tail(
            ::typeof(StatsFuns.$f),
            ::Val{T},
            k::Real,
            θ::Real,
            x::Real,
        ) where {T}
            primal = StatsFuns.$f(k, θ, x)
            logdensity = StatsFuns.gammalogpdf(k, θ, x)
            partial_θ, partial_x = _tail_partials(
                primal,
                logdensity,
                θ,
                x,
                $complementary,
                $logarithmic,
            )
            return primal, partial_θ, partial_x
        end

        ForwardDiff.@define_ternary_dual_op(
            StatsFuns.$f,
            _shape_is_not_differentiable(),
            _shape_is_not_differentiable(),
            _shape_is_not_differentiable(),
            begin
                vy, vz = ForwardDiff.value(Tyz, y), ForwardDiff.value(Tyz, z)
                primal, partial_θ, partial_x = _dual_tail(StatsFuns.$f, Val(Tyz), x, vy, vz)
                ForwardDiff.dual_definition_retval(
                    Val(Tyz),
                    primal,
                    partial_θ,
                    ForwardDiff.partials(Tyz, y),
                    partial_x,
                    ForwardDiff.partials(Tyz, z),
                )
            end,
            _shape_is_not_differentiable(),
            begin
                vy = ForwardDiff.value(Ty, y)
                primal, partial_θ, _ = _dual_tail(StatsFuns.$f, Val(Ty), x, vy, z)
                ForwardDiff.dual_definition_retval(
                    Val(Ty),
                    primal,
                    partial_θ,
                    ForwardDiff.partials(Ty, y),
                )
            end,
            begin
                vz = ForwardDiff.value(Tz, z)
                primal, _, partial_x = _dual_tail(StatsFuns.$f, Val(Tz), x, y, vz)
                ForwardDiff.dual_definition_retval(
                    Val(Tz),
                    primal,
                    partial_x,
                    ForwardDiff.partials(Tz, z),
                )
            end,
        )
    end
end

end # module
