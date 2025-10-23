# functions related to noncentral F distribution

# Rmath implementations
function nfdistpdf(k1::Real, k2::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k1, k2, λ, x))
    return convert(T, Rmath.dnf(x, k1, k2, λ, false))
end
function nfdistlogpdf(k1::Real, k2::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k1, k2, λ, x))
    return convert(T, Rmath.dnf(x, k1, k2, λ, true))
end

function nfdistcdf(k1::Real, k2::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k1, k2, λ, x))
    return convert(T, Rmath.pnf(x, k1, k2, λ, true, false))
end
function nfdistccdf(k1::Real, k2::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k1, k2, λ, x))
    return convert(T, Rmath.pnf(x, k1, k2, λ, false, false))
end
function nfdistlogcdf(k1::Real, k2::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k1, k2, λ, x))
    return convert(T, Rmath.pnf(x, k1, k2, λ, true, true))
end
function nfdistlogccdf(k1::Real, k2::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k1, k2, λ, x))
    return convert(T, Rmath.pnf(x, k1, k2, λ, false, true))
end

function nfdistinvcdf(k1::Real, k2::Real, λ::Real, q::Real)
    T = float(Base.promote_typeof(k1, k2, λ, q))
    return convert(T, Rmath.qnf(q, k1, k2, λ, true, false))
end
function nfdistinvccdf(k1::Real, k2::Real, λ::Real, q::Real)
    T = float(Base.promote_typeof(k1, k2, λ, q))
    return convert(T, Rmath.qnf(q, k1, k2, λ, false, false))
end
function nfdistinvlogcdf(k1::Real, k2::Real, λ::Real, lq::Real)
    T = float(Base.promote_typeof(k1, k2, λ, lq))
    return convert(T, Rmath.qnf(lq, k1, k2, λ, true, true))
end
function nfdistinvlogccdf(k1::Real, k2::Real, λ::Real, lq::Real)
    T = float(Base.promote_typeof(k1, k2, λ, lq))
    return convert(T, Rmath.qnf(lq, k1, k2, λ, false, true))
end
