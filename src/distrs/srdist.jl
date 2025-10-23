# functions related to studesrdistized range distribution

# Rmath implementations
function srdistcdf(k::Real, ν::Real, x::Real)
    T = float(Base.promote_typeof(k, ν, x))
    return convert(T, Rmath.ptukey(x, k, ν, true, false))
end
function srdistccdf(k::Real, ν::Real, x::Real)
    T = float(Base.promote_typeof(k, ν, x))
    return convert(T, Rmath.ptukey(x, k, ν, false, false))
end
function srdistlogcdf(k::Real, ν::Real, x::Real)
    T = float(Base.promote_typeof(k, ν, x))
    return convert(T, Rmath.ptukey(x, k, ν, true, true))
end
function srdistlogccdf(k::Real, ν::Real, x::Real)
    T = float(Base.promote_typeof(k, ν, x))
    return convert(T, Rmath.ptukey(x, k, ν, false, true))
end

function srdistinvcdf(k::Real, ν::Real, q::Real)
    T = float(Base.promote_typeof(k, ν, q))
    return convert(T, Rmath.qtukey(q, k, ν, true, false))
end
function srdistinvccdf(k::Real, ν::Real, q::Real)
    T = float(Base.promote_typeof(k, ν, q))
    return convert(T, Rmath.qtukey(q, k, ν, false, false))
end
function srdistinvlogcdf(k::Real, ν::Real, lq::Real)
    T = float(Base.promote_typeof(k, ν, lq))
    return convert(T, Rmath.qtukey(lq, k, ν, true, true))
end
function srdistinvlogccdf(k::Real, ν::Real, lq::Real)
    T = float(Base.promote_typeof(k, ν, lq))
    return convert(T, Rmath.qtukey(lq, k, ν, false, true))
end
