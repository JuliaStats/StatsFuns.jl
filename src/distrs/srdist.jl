# functions related to studentized range distribution

# Rmath implementations
# The `1` argument is `nranges` (number of ranges) which `Rmath.ptukey`/`qtukey`
# default to `1.0` in `kwargs` but is positional in the Julia wrapper. Without
# it the `lower_tail`/`log_p` Bool flags shift into the `nranges` slot.
function srdistcdf(k::Real, ν::Real, x::Real)
    T = float(Base.promote_typeof(k, ν, x))
    return convert(T, Rmath.ptukey(x, k, ν, 1, true, false))
end
function srdistccdf(k::Real, ν::Real, x::Real)
    T = float(Base.promote_typeof(k, ν, x))
    return convert(T, Rmath.ptukey(x, k, ν, 1, false, false))
end
function srdistlogcdf(k::Real, ν::Real, x::Real)
    T = float(Base.promote_typeof(k, ν, x))
    return convert(T, Rmath.ptukey(x, k, ν, 1, true, true))
end
function srdistlogccdf(k::Real, ν::Real, x::Real)
    T = float(Base.promote_typeof(k, ν, x))
    return convert(T, Rmath.ptukey(x, k, ν, 1, false, true))
end

function srdistinvcdf(k::Real, ν::Real, q::Real)
    T = float(Base.promote_typeof(k, ν, q))
    return convert(T, Rmath.qtukey(q, k, ν, 1, true, false))
end
function srdistinvccdf(k::Real, ν::Real, q::Real)
    T = float(Base.promote_typeof(k, ν, q))
    return convert(T, Rmath.qtukey(q, k, ν, 1, false, false))
end
function srdistinvlogcdf(k::Real, ν::Real, lq::Real)
    T = float(Base.promote_typeof(k, ν, lq))
    return convert(T, Rmath.qtukey(lq, k, ν, 1, true, true))
end
function srdistinvlogccdf(k::Real, ν::Real, lq::Real)
    T = float(Base.promote_typeof(k, ν, lq))
    return convert(T, Rmath.qtukey(lq, k, ν, 1, false, true))
end
