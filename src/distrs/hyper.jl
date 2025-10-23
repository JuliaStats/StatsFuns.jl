# functions related to hyper-geometric distribution

# Rmath implementations
function hyperpdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    return convert(T, Rmath.dhyper(x, ms, mf, n, false))
end
function hyperlogpdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    return convert(T, Rmath.dhyper(x, ms, mf, n, true))
end

function hypercdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    return convert(T, Rmath.phyper(x, ms, mf, n, true, false))
end
function hyperccdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    return convert(T, Rmath.phyper(x, ms, mf, n, false, false))
end
function hyperlogcdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    return convert(T, Rmath.phyper(x, ms, mf, n, true, true))
end
function hyperlogccdf(ms::Real, mf::Real, n::Real, x::Real)
    T = float(Base.promote_typeof(ms, mf, n, x))
    return convert(T, Rmath.phyper(x, ms, mf, n, false, true))
end

function hyperinvcdf(ms::Real, mf::Real, n::Real, q::Real)
    T = float(Base.promote_typeof(ms, mf, n, q))
    return convert(T, Rmath.qhyper(q, ms, mf, n, true, false))
end
function hyperinvccdf(ms::Real, mf::Real, n::Real, q::Real)
    T = float(Base.promote_typeof(ms, mf, n, q))
    return convert(T, Rmath.qhyper(q, ms, mf, n, false, false))
end
function hyperinvlogcdf(ms::Real, mf::Real, n::Real, lq::Real)
    T = float(Base.promote_typeof(ms, mf, n, lq))
    return convert(T, Rmath.qhyper(lq, ms, mf, n, true, true))
end
function hyperinvlogccdf(ms::Real, mf::Real, n::Real, lq::Real)
    T = float(Base.promote_typeof(ms, mf, n, lq))
    return convert(T, Rmath.qhyper(lq, ms, mf, n, false, true))
end
