# functions related to Poisson distribution

# Julia implementations
poispdf(λ::Real, x::Real) = exp(poislogpdf(λ, x))

function poislogpdf(λ::Real, x::Real)
    logupdf = poislogupdf(λ, x)
    return isfinite(logupdf) ? logupdf - λ : logupdf
end

poislogupdf(λ::Real, x::Real) = poislogupdf(promote(λ, x)...)
function poislogupdf(λ::T, x::T) where {T <: Real}
    val = xlogy(x, λ) - loggamma(x + 1)
    return x >= 0 && isinteger(x) ? val : oftype(val, -Inf)
end

function poislogulikelihood(λ::Real, x::Real)
    val = xlogy(x, λ) - λ
    return x >= 0 && isinteger(x) ? val : oftype(val, -Inf)
end

# Just use the Gamma definitions
poiscdf(λ::Real, x::Real) = gammaccdf(max(0, floor(x + 1)), 1, λ)

poisccdf(λ::Real, x::Real) = gammacdf(max(0, floor(x + 1)), 1, λ)

poislogcdf(λ::Real, x::Real) = gammalogccdf(max(0, floor(x + 1)), 1, λ)

poislogccdf(λ::Real, x::Real) = gammalogcdf(max(0, floor(x + 1)), 1, λ)

# Rmath implementations
function poisinvcdf(λ::Real, q::Real)
    T = float(Base.promote_typeof(λ, q))
    return convert(T, Rmath.qpois(q, λ, true, false))
end
function poisinvccdf(λ::Real, q::Real)
    T = float(Base.promote_typeof(λ, q))
    return convert(T, Rmath.qpois(q, λ, false, false))
end
function poisinvlogcdf(λ::Real, lq::Real)
    T = float(Base.promote_typeof(λ, lq))
    return convert(T, Rmath.qpois(lq, λ, true, true))
end
function poisinvlogccdf(λ::Real, lq::Real)
    T = float(Base.promote_typeof(λ, lq))
    return convert(T, Rmath.qpois(lq, λ, false, true))
end
