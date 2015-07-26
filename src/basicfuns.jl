# common facilities

f64(x::Real) = convert(Float64, x)

# scalar functions

xlogx(x::FloatingPoint) = x > zero(x) ? x * log(x) : zero(x)
xlogx(x::Real) = xlogx(float(x))

xlogy{T<:FloatingPoint}(x::T, y::T) = x > zero(T) ? x * log(y) : zero(x)
xlogy{T<:Real}(x::T, y::T) = xlogy(float(x), float(y))
xlogy(x::Real, y::Real) = xlogy(promote(x, y)...)

# logistic: 1 / (1 + exp(-x))
#
logistic(x::FloatingPoint) = one(x) / (one(x) + exp(-x))
logistic(x::Real) = logistic(float(x))

# logit: log(x / (1 - x))
#
logit(x::FloatingPoint) = log(x / (one(x) - x))
logit(x::Real) = logit(float(x))

# log1psq: log(1+x^2)
#
log1psq(x::FloatingPoint) = log1p(abs2(x))
log1psq(x::Union(Float32,Float64)) = (ax = abs(x); ax < maxintfloat(x) ? log1p(abs2(ax)) : 2 * log(ax))
log1psq(x::Real) = log1psq(float(x))

# log1pexp: log(1+exp(x))
#
log1pexp(x::FloatingPoint) = x < 18.0 ? log1p(exp(x)) : x < 33.3 ? x + exp(-x) : x
log1pexp(x::Float32) = x < 9.0f0 ? log1p(exp(x)) : x < 16.0f0 ? x + exp(-x) : x
log1pexp(x::Real) = log1pexp(float(x))

# log1mexp: log(1 - exp(x))
#
# See:
#   Martin Maechler (2012) "Accurately Computing log(1 − exp(− |a|))"
#   http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
#
# Note: different than Maechler (2012), no negation inside parantheses
log1mexp(x::FloatingPoint) = x < loghalf ? log1p(-exp(x)) : log(-expm1(x))

# log2mexp: log(2 - exp(x))
#
log2mexp(x::FloatingPoint) = log1p(-expm1(x))
log2mexp(x::Real) = log2mexp(float(x))

# logexpm1: log(exp(x) - 1)
#
logexpm1(x::FloatingPoint) = x <= 18.0 ? log(expm1(x)) : x <= 33.3 ? x - exp(-x) : x
logexpm1(x::Float32) = x <= 9f0 ? log(expm1(x)) : x <= 16f0 ? x - exp(-x) : x
logexpm1(x::Real) = logexpm1(float(x))

@vectorize_1arg Real xlogx
@vectorize_2arg Real xlogy
@vectorize_1arg Real logistic
@vectorize_1arg Real logit
@vectorize_1arg Real log1psq
@vectorize_1arg Real log1pexp
@vectorize_1arg Real log1mexp
@vectorize_1arg Real log2mexp
@vectorize_1arg Real logexpm1

const softplus = log1pexp
const invsoftplus = logexpm1

# log1pmx: log(1 + x) - x
#
# use naive calculation or range reduction outside kernel range.
# accurate ~2ulps for all x
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

# logmxp1: log(x) - x + 1
#
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


## logsumexp

logsumexp{T<:FloatingPoint}(x::T, y::T) = x > y ? x + log1p(exp(y - x)) : y + log1p(exp(x - y))
logsumexp{T<:Real}(x::T, y::T) = logsumexp(float(x), float(y))
logsumexp(x::Real, y::Real) = logsumexp(promote(x, y)...)

function logsumexp{T<:Real}(x::AbstractArray{T})
    isempty(x) && return -Inf
    u = maximum(x)
    s = 0.
    for i = 1:length(x)
    	@inbounds s += exp(x[i] - u)
    end
    log(s) + u
end

## softmax

function softmax!{R<:FloatingPoint,T<:Real}(r::AbstractArray{R}, x::AbstractArray{T})
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

softmax!{T<:FloatingPoint}(x::AbstractArray{T}) = softmax!(x, x)
softmax{T<:Real}(x::AbstractArray{T}) = softmax!(Array(Float64, size(x)), x)
