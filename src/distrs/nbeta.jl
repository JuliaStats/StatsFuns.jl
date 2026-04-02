# functions related to noncentral beta distribution

# Rmath implementations
function nbetapdf(α::Real, β::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(α, β, λ, x))
    return convert(T, Rmath.dnbeta(x, α, β, λ, false))
end
function nbetalogpdf(α::Real, β::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(α, β, λ, x))
    return convert(T, Rmath.dnbeta(x, α, β, λ, true))
end

function nbetacdf(α::Real, β::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(α, β, λ, x))
    return convert(T, Rmath.pnbeta(x, α, β, λ, true, false))
end
function nbetaccdf(α::Real, β::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(α, β, λ, x))
    return convert(T, Rmath.pnbeta(x, α, β, λ, false, false))
end
function nbetalogcdf(α::Real, β::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(α, β, λ, x))
    return convert(T, Rmath.pnbeta(x, α, β, λ, true, true))
end
function nbetalogccdf(α::Real, β::Real, λ::Real, x::Real)
    T = float(Base.promote_typeof(α, β, λ, x))
    return convert(T, Rmath.pnbeta(x, α, β, λ, false, true))
end

function nbetainvcdf(α::Real, β::Real, λ::Real, q::Real)
    T = float(Base.promote_typeof(α, β, λ, q))
    return convert(T, Rmath.qnbeta(q, α, β, λ, true, false))
end
function nbetainvccdf(α::Real, β::Real, λ::Real, q::Real)
    T = float(Base.promote_typeof(α, β, λ, q))
    return convert(T, Rmath.qnbeta(q, α, β, λ, false, false))
end
function nbetainvlogcdf(α::Real, β::Real, λ::Real, lq::Real)
    T = float(Base.promote_typeof(α, β, λ, lq))
    return convert(T, Rmath.qnbeta(lq, α, β, λ, true, true))
end
function nbetainvlogccdf(α::Real, β::Real, λ::Real, lq::Real)
    T = float(Base.promote_typeof(α, β, λ, lq))
    return convert(T, Rmath.qnbeta(lq, α, β, λ, false, true))
end
