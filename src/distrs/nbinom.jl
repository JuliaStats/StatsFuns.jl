# functions related to negative binomial distribution

# Rmath implementations
function nbinompdf(r::Real, p::Real, x::Real)
    T = float(Base.promote_typeof(r, p, x))
    return convert(T, Rmath.dnbinom(x, r, p, false))
end
function nbinomlogpdf(r::Real, p::Real, x::Real)
    T = float(Base.promote_typeof(r, p, x))
    return convert(T, Rmath.dnbinom(x, r, p, true))
end

function nbinomcdf(r::Real, p::Real, x::Real)
    T = float(Base.promote_typeof(r, p, x))
    return convert(T, Rmath.pnbinom(x, r, p, true, false))
end
function nbinomccdf(r::Real, p::Real, x::Real)
    T = float(Base.promote_typeof(r, p, x))
    return convert(T, Rmath.pnbinom(x, r, p, false, false))
end
function nbinomlogcdf(r::Real, p::Real, x::Real)
    T = float(Base.promote_typeof(r, p, x))
    return convert(T, Rmath.pnbinom(x, r, p, true, true))
end
function nbinomlogccdf(r::Real, p::Real, x::Real)
    T = float(Base.promote_typeof(r, p, x))
    return convert(T, Rmath.pnbinom(x, r, p, false, true))
end

function nbinominvcdf(r::Real, p::Real, q::Real)
    T = float(Base.promote_typeof(r, p, q))
    return convert(T, Rmath.qnbinom(q, r, p, true, false))
end
function nbinominvccdf(r::Real, p::Real, q::Real)
    T = float(Base.promote_typeof(r, p, q))
    return convert(T, Rmath.qnbinom(q, r, p, false, false))
end
function nbinominvlogcdf(r::Real, p::Real, lq::Real)
    T = float(Base.promote_typeof(r, p, lq))
    return convert(T, Rmath.qnbinom(lq, r, p, true, true))
end
function nbinominvlogccdf(r::Real, p::Real, lq::Real)
    T = float(Base.promote_typeof(r, p, lq))
    return convert(T, Rmath.qnbinom(lq, r, p, false, true))
end
