# functions related to Poisson distribution

# Julia implementations
poispdf(λ::Real, x::Real) = exp(poislogpdf(λ, x))

poislogpdf(λ::Real, x::Real) = poislogpdf(promote(λ, x)...)
function poislogpdf(λ::T, x::T) where {T <: Real}
    val = xlogy(x, λ) - λ - loggamma(x + 1)
    return x >= 0 && isinteger(x) ? val : oftype(val, -Inf)
end

# Just use the Gamma definitions
poiscdf(λ::Real, x::Real) = gammaccdf(max(0, floor(x + 1)), 1, λ)

poisccdf(λ::Real, x::Real) = gammacdf(max(0, floor(x + 1)), 1, λ)

poislogcdf(λ::Real, x::Real) = gammalogccdf(max(0, floor(x + 1)), 1, λ)

poislogccdf(λ::Real, x::Real) = gammalogcdf(max(0, floor(x + 1)), 1, λ)

# Inverse CDF: find smallest k such that poiscdf(λ, k) >= cprob.
# Based on VBA critpoiss by Ian Smith.
# PMF ratio: PMF(k+1)/PMF(k) = λ/(k+1), PMF(k-1)/PMF(k) = k/λ
function _critpoiss(λ::Float64, cprob::Float64)
    # Normal approximation for initial guess
    σ = sqrt(λ)
    i = max(0.0, floor(λ + norminvcdf(min(cprob, 1.0 - 1.0e-15)) * σ + 0.5))

    # Compute CDF at the guess
    pr = Float64(poiscdf(λ, i))

    if pr >= cprob
        # Search left: find smallest k with CDF(k) >= cprob
        while i > 0
            tpr = Float64(poispdf(λ, i))
            if pr - tpr < cprob
                return i
            end
            pr -= tpr
            i -= 1.0
        end
        return 0.0
    else
        # Search right (with guard against infinite loop)
        for _ in 1:10_000
            i += 1.0
            tpr = Float64(poispdf(λ, i))
            pr += tpr
            if pr >= cprob
                return i
            end
        end
        return i
    end
end

# Inverse CCDF: find smallest k such that poisccdf(λ, k) <= cprob.
# Based on VBA critcomppoiss by Ian Smith.
function _critcomppoiss(λ::Float64, cprob::Float64)
    # Normal approximation
    σ = sqrt(λ)
    i = max(0.0, floor(λ - norminvcdf(min(cprob, 1.0 - 1.0e-15)) * σ + 0.5))

    # Compute CCDF at the guess
    pr = Float64(poisccdf(λ, i))

    if pr > cprob
        # Search right (with guard)
        for _ in 1:10_000
            i += 1.0
            tpr = Float64(poispdf(λ, i))
            pr -= tpr
            if pr <= cprob
                return i
            end
        end
        return i
    else
        # Search left
        while i > 0
            tpr = Float64(poispdf(λ, i))
            if pr + tpr > cprob
                return i
            end
            pr += tpr
            i -= 1.0
        end
        return 0.0
    end
end

# Wrapper with edge cases
function _pois_invcdf(λ::Float64, q::Float64)
    if q < 0 || q > 1 || λ < 0 || isnan(q) || isnan(λ)
        return NaN
    elseif q == 0 || λ == 0
        return 0.0
    elseif q == 1
        return Inf
    end

    i = _critpoiss(λ, q)

    # Post-correction
    pr = Float64(poiscdf(λ, i))
    if pr >= q
        while i > 0
            pr2 = Float64(poiscdf(λ, i - 1.0))
            if pr2 < q
                return i
            end
            i -= 1.0
        end
        return 0.0
    else
        return i + 1.0
    end
end

function _pois_invccdf(λ::Float64, q::Float64)
    if q < 0 || q > 1 || λ < 0 || isnan(q) || isnan(λ)
        return NaN
    elseif q == 0
        return Inf
    elseif λ == 0
        return 0.0
    end

    i = _critcomppoiss(λ, q)

    # Post-correction
    pr = Float64(poisccdf(λ, i))
    if pr <= q
        while i > 0
            pr2 = Float64(poisccdf(λ, i - 1.0))
            if pr2 > q
                return i
            end
            i -= 1.0
        end
        return 0.0
    else
        return i + 1.0
    end
end

# Public API

function poisinvcdf(λ::Real, q::Real)
    T = float(Base.promote_typeof(λ, q))
    return convert(T, _pois_invcdf(Float64(λ), Float64(q)))
end

function poisinvccdf(λ::Real, q::Real)
    T = float(Base.promote_typeof(λ, q))
    return convert(T, _pois_invccdf(Float64(λ), Float64(q)))
end

function poisinvlogcdf(λ::Real, lq::Real)
    T = float(Base.promote_typeof(λ, lq))
    _lq = Float64(lq)
    # Use -expm1(lq) = 1-exp(lq) via invccdf for accuracy when lq ≈ 0
    result = if _lq > -1
        _pois_invccdf(Float64(λ), -expm1(_lq))
    else
        _pois_invcdf(Float64(λ), exp(_lq))
    end
    return convert(T, result)
end

function poisinvlogccdf(λ::Real, lq::Real)
    T = float(Base.promote_typeof(λ, lq))
    _lq = Float64(lq)
    result = if _lq > -1
        _pois_invcdf(Float64(λ), -expm1(_lq))
    else
        _pois_invccdf(Float64(λ), exp(_lq))
    end
    return convert(T, result)
end
