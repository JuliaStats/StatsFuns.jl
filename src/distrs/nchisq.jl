# functions related to noncentral chi-square distribution

# Rmath implementations
function nchisqpdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.dnchisq(x, k, λ, false))
end
function nchisqlogpdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.dnchisq(x, k, λ, true))
end

function nchisqcdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.pnchisq(x, k, λ, true, false))
end
function nchisqccdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.pnchisq(x, k, λ, false, false))
end
function nchisqlogcdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.pnchisq(x, k, λ, true, true))
end
function nchisqlogccdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.pnchisq(x, k, λ, false, true))
end

function nchisqinvcdf(k::Real, λ::Real, q::Real)
    T = float(Base.promote_typeof(k, λ, q))
    return convert(T, Rmath.qnchisq(q, k, λ, true, false))
end
function nchisqinvccdf(k::Real, λ::Real, q::Real)
    T = float(Base.promote_typeof(k, λ, q))
    return convert(T, Rmath.qnchisq(q, k, λ, false, false))
end
function nchisqinvlogcdf(k::Real, λ::Real, lq::Real)
    T = float(Base.promote_typeof(k, λ, lq))
    return convert(T, Rmath.qnchisq(lq, k, λ, true, true))
end
function nchisqinvlogccdf(k::Real, λ::Real, lq::Real)
    T = float(Base.promote_typeof(k, λ, lq))
    return convert(T, Rmath.qnchisq(lq, k, λ, false, true))
end
