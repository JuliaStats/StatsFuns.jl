module ForwardDiffExt

import StatsFuns
using ForwardDiff: ForwardDiff, ≺ # `≺` is used by `ForwardDiff.@define_ternary_dual_op`

# `SpecialFunctions.gamma_inc`, which every gamma tail bottoms out in, is defined
# for `Float16`, `Float32` and `Float64` only, so for `Dual` arguments we evaluate
# the function on the primal values and attach the partials analytically.

# With the tail written as `P(k, x/θ)`, `|∂P/∂x| = f_k(x)` and `|∂P/∂θ| = k f_{k+1}(x)`,
# both in log space so they are zero outside the support and at `x = 0`. The log forms
# divide them by the tail probability `exp(primal)`.
function _dual_tail(f, complementary::Bool, logarithmic::Bool, k::Real, θ::Real, x::Real)
    primal = f(k, θ, x)
    log_partial_x = StatsFuns.gammalogpdf(k, θ, x)
    log_partial_θ = log(k) + StatsFuns.gammalogpdf(k + 1, θ, x)
    sign = complementary ? -1 : 1
    partial_θ = -sign * _partial_magnitude(log_partial_θ, primal, logarithmic)
    partial_x = sign * _partial_magnitude(log_partial_x, primal, logarithmic)
    return primal, partial_θ, partial_x
end

function _partial_magnitude(log_partial::Real, primal::Real, logarithmic::Bool)
    if log_partial == -Inf
        return zero(log_partial)
    end
    return exp(logarithmic ? log_partial - primal : log_partial)
end

# `∂/∂k` of the regularized incomplete gamma has no elementary closed form. A `Dual`
# shape with zero partials, e.g. one promoted along with `θ` by `Distributions.Gamma`,
# is replaced by its value.
function _shape_value(::Type{T}, k::Real) where {T}
    if !iszero(ForwardDiff.partials(T, k))
        throw(
            ArgumentError(
                "the gamma tail functions are not differentiable with respect to " *
                    "the shape parameter",
            )
        )
    end
    return ForwardDiff.value(T, k)
end

for (f, complementary, logarithmic) in (
        (:gammacdf, false, false),
        (:gammaccdf, true, false),
        (:gammalogcdf, false, true),
        (:gammalogccdf, true, true),
    )
    @eval ForwardDiff.@define_ternary_dual_op(
        StatsFuns.$f,
        StatsFuns.$f(_shape_value(Txyz, x), y, z),
        StatsFuns.$f(_shape_value(Txy, x), y, z),
        StatsFuns.$f(_shape_value(Txz, x), y, z),
        begin
            vy, vz = ForwardDiff.value(Tyz, y), ForwardDiff.value(Tyz, z)
            primal, partial_θ, partial_x = _dual_tail(
                StatsFuns.$f,
                $complementary,
                $logarithmic,
                x,
                vy,
                vz,
            )
            ForwardDiff.dual_definition_retval(
                Val(Tyz),
                primal,
                partial_θ,
                ForwardDiff.partials(Tyz, y),
                partial_x,
                ForwardDiff.partials(Tyz, z),
            )
        end,
        StatsFuns.$f(_shape_value(Tx, x), y, z),
        begin
            vy = ForwardDiff.value(Ty, y)
            primal, partial_θ, _ = _dual_tail(
                StatsFuns.$f,
                $complementary,
                $logarithmic,
                x,
                vy,
                z,
            )
            ForwardDiff.dual_definition_retval(
                Val(Ty),
                primal,
                partial_θ,
                ForwardDiff.partials(Ty, y),
            )
        end,
        begin
            vz = ForwardDiff.value(Tz, z)
            primal, _, partial_x = _dual_tail(
                StatsFuns.$f,
                $complementary,
                $logarithmic,
                x,
                y,
                vz,
            )
            ForwardDiff.dual_definition_retval(
                Val(Tz),
                primal,
                partial_x,
                ForwardDiff.partials(Tz, z),
            )
        end,
    )
end

end # module
