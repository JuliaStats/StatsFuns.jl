# functions related to noncentral T distribution

ntdistlogupdf(k::Real, λ::Real, x::Real) = ntdistlogpdf(k, λ, x)
ntdistlogulikelihood(k::Real, λ::Real, x::Real) = ntdistlogpdf(k, λ, x)

# Rmath implementations
function ntdistpdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.dnt(x, k, λ, false))
end
function ntdistlogpdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.dnt(x, k, λ, true))
end

function ntdistcdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.pnt(x, k, λ, true, false))
end
function ntdistccdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.pnt(x, k, λ, false, false))
end
function ntdistlogcdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.pnt(x, k, λ, true, true))
end
function ntdistlogccdf(k::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(k, λ, x))
    return convert(T, Rmath.pnt(x, k, λ, false, true))
end

function ntdistinvcdf(k::Real, λ::Real, q::Real)
    T = float(Base.promote_typeof(k, λ, q))
    return convert(T, Rmath.qnt(q, k, λ, true, false))
end
function ntdistinvccdf(k::Real, λ::Real, q::Real)
    T = float(Base.promote_typeof(k, λ, q))
    return convert(T, Rmath.qnt(q, k, λ, false, false))
end
function ntdistinvlogcdf(k::Real, λ::Real, lq::Real)
    T = float(Base.promote_typeof(k, λ, lq))
    return convert(T, Rmath.qnt(lq, k, λ, true, true))
end
function ntdistinvlogccdf(k::Real, λ::Real, lq::Real)
    T = float(Base.promote_typeof(k, λ, lq))
    return convert(T, Rmath.qnt(lq, k, λ, false, true))
end
