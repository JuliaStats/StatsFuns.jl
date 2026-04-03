# functions related to binomial distribution

# Julia implementations
binompdf(n::Real, p::Real, k::Real) = exp(binomlogpdf(n, p, k))

binomlogpdf(n::Real, p::Real, k::Real) = binomlogpdf(promote(n, p, k)...)
function binomlogpdf(n::T, p::T, k::T) where {T <: Real}
    m = clamp(k, 0, n)
    val = min(0, betalogpdf(m + 1, n - m + 1, p) - log(n + 1))
    return 0 <= k <= n && isinteger(k) ? val : oftype(val, -Inf)
end

for l in ("", "log"), compl in (false, true)
    fbinom = Symbol(string("binom", l, ifelse(compl, "c", ""), "cdf"))
    fbeta = Symbol(string("beta", l, ifelse(compl, "", "c"), "cdf"))
    @eval function ($fbinom)(n::Real, p::Real, k::Real)
        if isnan(k)
            return last(promote(n, p, k))
        end
        res = ($fbeta)(max(0, floor(k) + 1), max(0, n - floor(k)), p)

        # When p == 1 the betaccdf doesn't return the correct result
        # so these cases have to be special cased
        if isone(p)
            newres = oftype(res, $compl ? k < n : k >= n)
            return $(l === "" ? :newres : :(log(newres)))
        else
            return res
        end
    end
end

# Inverse CDF: find smallest k such that binomcdf(n, p, k) >= cprob.
# Based on VBA critbinomial by Ian Smith.
# PMF ratios: PMF(k+1)/PMF(k) = (n-k)*p/((k+1)*(1-p))
#             PMF(k-1)/PMF(k) = k*(1-p)/((n-k+1)*p)
function _critbinomial(n::Float64, p::Float64, cprob::Float64)
    q = 1.0 - p

    # Normal approximation for initial guess
    σ = sqrt(n * p * q)
    i = clamp(floor(n * p + norminvcdf(min(cprob, 1.0 - 1.0e-15)) * σ + 0.5), 0.0, n)

    # Compute CDF at the guess
    pr = Float64(binomcdf(n, p, i))

    if pr >= cprob
        # Search left: find smallest k with CDF(k) >= cprob
        while i > 0
            tpr = Float64(binompdf(n, p, i))
            if pr - tpr < cprob
                return i  # CDF(i-1) < cprob <= CDF(i)
            end
            pr -= tpr
            i -= 1.0
        end
        return 0.0
    else
        # Search right
        while i < n
            i += 1.0
            tpr = Float64(binompdf(n, p, i))
            pr += tpr
            if pr >= cprob
                return i
            end
        end
        return n
    end
end

# Inverse CCDF: find smallest k such that binomccdf(n, p, k) <= cprob.
# Based on VBA critcompbinomial by Ian Smith.
function _critcompbinomial(n::Float64, p::Float64, cprob::Float64)
    q = 1.0 - p

    # Normal approximation
    σ = sqrt(n * p * q)
    i = clamp(floor(n * p - norminvcdf(min(cprob, 1.0 - 1.0e-15)) * σ + 0.5), 0.0, n)

    # Compute CCDF at the guess
    pr = Float64(binomccdf(n, p, i))

    if pr > cprob
        # Search right
        while i < n
            i += 1.0
            tpr = Float64(binompdf(n, p, i))
            pr -= tpr
            if pr <= cprob
                return i
            end
        end
        return n
    else
        # Search left: find smallest k with CCDF(k) <= cprob
        while i > 0
            tpr = Float64(binompdf(n, p, i))
            if pr + tpr > cprob
                return i  # CCDF(i) <= cprob < CCDF(i-1)
            end
            pr += tpr
            i -= 1.0
        end
        return 0.0
    end
end

# Wrapper with edge cases and post-correction (VBA crit_binomial)
function _binom_invcdf(n::Float64, p::Float64, q::Float64)
    if q < 0 || q > 1 || n < 0 || p < 0 || p > 1 || isnan(q) || isnan(n) || isnan(p)
        return NaN
    elseif q == 0 || p == 0
        return 0.0
    elseif p == 1
        return n
    end

    i = _critbinomial(n, p, q)

    # Post-correction
    pr = Float64(binomcdf(n, p, i))
    if pr >= q
        while i > 0
            pr2 = Float64(binomcdf(n, p, i - 1.0))
            if pr2 < q
                return i
            end
            i -= 1.0
            pr = pr2
        end
        return 0.0
    else
        return i + 1.0
    end
end

# Wrapper with edge cases and post-correction (VBA comp_crit_binomial)
function _binom_invccdf(n::Float64, p::Float64, q::Float64)
    if q < 0 || q > 1 || n < 0 || p < 0 || p > 1 || isnan(q) || isnan(n) || isnan(p)
        return NaN
    elseif q == 1 || p == 0
        return 0.0
    elseif q == 0 || p == 1
        return n
    end

    i = _critcompbinomial(n, p, q)

    # Post-correction
    pr = Float64(binomccdf(n, p, i))
    if pr <= q
        while i > 0
            pr2 = Float64(binomccdf(n, p, i - 1.0))
            if pr2 > q
                return i
            end
            i -= 1.0
            pr = pr2
        end
        return 0.0
    else
        return i + 1.0
    end
end

# Public API

function binominvcdf(n::Real, p::Real, q::Real)
    T = float(Base.promote_typeof(n, p, q))
    return convert(T, _binom_invcdf(Float64(n), Float64(p), Float64(q)))
end

function binominvccdf(n::Real, p::Real, q::Real)
    T = float(Base.promote_typeof(n, p, q))
    return convert(T, _binom_invccdf(Float64(n), Float64(p), Float64(q)))
end

function binominvlogcdf(n::Real, p::Real, lq::Real)
    T = float(Base.promote_typeof(n, p, lq))
    return convert(T, _binom_invcdf(Float64(n), Float64(p), exp(Float64(lq))))
end

function binominvlogccdf(n::Real, p::Real, lq::Real)
    T = float(Base.promote_typeof(n, p, lq))
    return convert(T, _binom_invccdf(Float64(n), Float64(p), exp(Float64(lq))))
end
