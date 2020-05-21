# common facilities

# scalar functions
"""
    xlogx(x::Real)

Compute `x * log(x)`, returning zero if `x = 0`.

```jldoctest
julia> StatsFuns.xlogx(0)
0.0
```
"""
function xlogx(x::Number)
    result = x * log(x)
    ifelse(iszero(x), zero(result), result)
end

"""
    xlogy(x::Real, y::Real)

Compute `x * log(y)`, returning zero if `x = 0`.

```jldoctest
julia> StatsFuns.xlogy(0, 0)
0.0
```
"""
function xlogy(x::Number, y::Number)
    result = x * log(y)
    ifelse(iszero(x) && !isnan(y), zero(result), result)
end

# The following bounds are precomputed versions of the following abstract
# function, but the implicit interface for AbstractFloat doesn't uniformly
# enforce that all floating point types implement nextfloat and prevfloat.
# @inline function _logistic_bounds(x::AbstractFloat)
#     (
#         logit(nextfloat(zero(float(x)))),
#         logit(prevfloat(one(float(x)))),
#     )
# end

@inline _logistic_bounds(x::Float16) = (Float16(-16.64), Float16(7.625))
@inline _logistic_bounds(x::Float32) = (-103.27893f0, 16.635532f0)
@inline _logistic_bounds(x::Float64) = (-744.4400719213812, 36.7368005696771)

"""
    logistic(x::Real)

The [logistic](https://en.wikipedia.org/wiki/Logistic_function) sigmoid function mapping a real number to a value in the interval [0,1],
```math
\\sigma(x) = \\frac{1}{e^{-x} + 1} = \\frac{e^x}{1+e^x}.
```

Its inverse is the [`logit`](@ref) function.
"""
logistic(x::Real) = inv(exp(-x) + one(x))

function logistic(x::Union{Float16, Float32, Float64})
    e = exp(x)
    lower, upper = _logistic_bounds(x)
    ifelse(
        x < lower,
        zero(x),
        ifelse(
            x > upper,
            one(x),
            e / (one(x) + e)
        )
    )
end

"""
    logit(p::Real)

The [logit](https://en.wikipedia.org/wiki/Logit) or log-odds transformation,
```math
\\log\\left(\\frac{x}{1-x}\\right), \\text{where} 0 < x < 1
```
Its inverse is the [`logistic`](@ref) function.
"""
logit(x::Real) = log(x / (one(x) - x))

"""
    log1psq(x::Real)

Return `log(1+x^2)` evaluated carefully for `abs(x)` very small or very large.
"""
log1psq(x::Real) = log1p(abs2(x))
function log1psq(x::Union{Float32,Float64})
    ax = abs(x)
    ax < maxintfloat(x) ? log1p(abs2(ax)) : 2 * log(ax)
end

"""
    log1pexp(x::Real)

Return `log(1+exp(x))` evaluated carefully for largish `x`.

This is also called the ["softplus"](https://en.wikipedia.org/wiki/Rectifier_(neural_networks))
transformation, being a smooth approximation to `max(0,x)`. Its inverse is [`logexpm1`](@ref).
"""
log1pexp(x::Real) = x < 18.0 ? log1p(exp(x)) : x < 33.3 ? x + exp(-x) : oftype(exp(-x), x)
log1pexp(x::Float32) = x < 9.0f0 ? log1p(exp(x)) : x < 16.0f0 ? x + exp(-x) : oftype(exp(-x), x)

"""
    log1mexp(x::Real)

Return `log(1 - exp(x))`

See:
 * Martin Maechler (2012) "Accurately Computing log(1 − exp(− |a|))",
   http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf

Note: different than Maechler (2012), no negation inside parentheses
"""
log1mexp(x::Real) = x < loghalf ? log1p(-exp(x)) : log(-expm1(x))

"""
    log2mexp(x::Real)

Return `log(2 - exp(x))` evaluated as `log1p(-expm1(x))`
"""
log2mexp(x::Real) = log1p(-expm1(x))

"""
    logexpm1(x::Real)

Return `log(exp(x) - 1)` or the "invsoftplus" function.
It is the inverse of [`log1pexp`](@ref) (aka "softplus").
"""
logexpm1(x::Real) = x <= 18.0 ? log(expm1(x)) : x <= 33.3 ? x - exp(-x) : oftype(exp(-x), x)
logexpm1(x::Float32) = x <= 9f0 ? log(expm1(x)) : x <= 16f0 ? x - exp(-x) : oftype(exp(-x), x)

const softplus = log1pexp
const invsoftplus = logexpm1

"""
    log1pmx(x::Float64)

Return `log(1 + x) - x`

Use naive calculation or range reduction outside kernel range.  Accurate ~2ulps for all `x`.
"""
function log1pmx(x::Float64)
    if !(-0.7 < x < 0.9)
        return log1p(x) - x
    elseif x > 0.315
        u = (x-0.5)/1.5
        return _log1pmx_ker(u) - 9.45348918918356180e-2 - 0.5*u
    elseif x > -0.227
        return _log1pmx_ker(x)
    elseif x > -0.4
        u = (x+0.25)/0.75
        return _log1pmx_ker(u) - 3.76820724517809274e-2 + 0.25*u
    elseif x > -0.6
        u = (x+0.5)*2.0
        return _log1pmx_ker(u) - 1.93147180559945309e-1 + 0.5*u
    else
        u = (x+0.625)/0.375
        return _log1pmx_ker(u) - 3.55829253011726237e-1 + 0.625*u
    end
end

"""
    logmxp1(x::Float64)

Return `log(x) - x + 1` carefully evaluated.
"""
function logmxp1(x::Float64)
    if x <= 0.3
        return (log(x) + 1.0) - x
    elseif x <= 0.4
        u = (x-0.375)/0.375
        return _log1pmx_ker(u) - 3.55829253011726237e-1 + 0.625*u
    elseif x <= 0.6
        u = 2.0*(x-0.5)
        return _log1pmx_ker(u) - 1.93147180559945309e-1 + 0.5*u
    else
        return log1pmx(x - 1.0)
    end
end

# The kernel of log1pmx
# Accuracy within ~2ulps for -0.227 < x < 0.315
function _log1pmx_ker(x::Float64)
    r = x/(x+2.0)
    t = r*r
    w = @horner(t,
                6.66666666666666667e-1, # 2/3
                4.00000000000000000e-1, # 2/5
                2.85714285714285714e-1, # 2/7
                2.22222222222222222e-1, # 2/9
                1.81818181818181818e-1, # 2/11
                1.53846153846153846e-1, # 2/13
                1.33333333333333333e-1, # 2/15
                1.17647058823529412e-1) # 2/17
    hxsq = 0.5*x*x
    r*(hxsq+w*t)-hxsq
end


"""
    logaddexp(x::Real, y::Real)

Return `log(exp(x) + exp(y))`, avoiding intermediate overflow/undeflow, and handling non-finite values.
"""
function logaddexp(x::Real, y::Real)
    # ensure Δ = 0 if x = y = Inf
    Δ = ifelse(x == y, zero(x - y), abs(x - y))
    max(x, y) + log1pexp(-Δ)
end

Base.@deprecate logsumexp(x::Real, y::Real) logaddexp(x, y)

"""
    logsubexp(x, y)

Return `log(abs(e^x - e^y))`, preserving numerical accuracy.
"""
logsubexp(x::Real, y::Real) = max(x, y) + log1mexp(-abs(x - y))

"""
    logsumexp(X)

Compute `log(sum(exp, X))`, evaluated avoiding intermediate overflow/undeflow.

`X` should be an iterator of real numbers.
"""
function logsumexp(X)
    isempty(X) && return log(sum(X))
    reduce(logaddexp, X)
end
function logsumexp(X::AbstractArray{T}; dims=:) where {T<:Real}
    # Do not use log(zero(T)) directly to avoid issues with ForwardDiff (#82)
    u = reduce(max, X, dims=dims, init=oftype(log(zero(T)), -Inf))
    u isa AbstractArray || isfinite(u) || return float(u)
    let u=u # avoid https://github.com/JuliaLang/julia/issues/15276
        # TODO: remove the branch when JuliaLang/julia#31020 is merged.
        if u isa AbstractArray
            u .+ log.(sum(exp.(X .- u); dims=dims))
        else
            u + log(sum(x -> exp(x-u), X))
        end
    end
end


"""
    softmax!(r::AbstractArray, x::AbstractArray)

Overwrite `r` with the `softmax` (or _normalized exponential_) transformation of `x`

That is, `r` is overwritten with `exp.(x)`, normalized to sum to 1.

See the [Wikipedia entry](https://en.wikipedia.org/wiki/Softmax_function)
"""
function softmax!(r::AbstractArray{R}, x::AbstractArray{T}) where {R<:AbstractFloat,T<:Real}
    n = length(x)
    length(r) == n || throw(DimensionMismatch("Inconsistent array lengths."))
    u = maximum(x)
    s = 0.
    @inbounds for i = 1:n
        s += (r[i] = exp(x[i] - u))
    end
    invs = convert(R, inv(s))
    @inbounds for i = 1:n
        r[i] *= invs
    end
    r
end

"""
    softmax(x::AbstractArray{<:Real})

Return the [`softmax transformation`](https://en.wikipedia.org/wiki/Softmax_function) applied to `x`
"""
softmax!(x::AbstractArray{<:AbstractFloat}) = softmax!(x, x)
softmax(x::AbstractArray{<:Real}) = softmax!(similar(x, Float64), x)
