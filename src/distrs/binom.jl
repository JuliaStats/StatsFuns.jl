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

# Rmath implementations
function binominvcdf(n::Real, p::Real, q::Real)
    T = float(Base.promote_typeof(n, p, q))
    return convert(T, Rmath.qbinom(q, n, p, true, false))
end
function binominvccdf(n::Real, p::Real, q::Real)
    T = float(Base.promote_typeof(n, p, q))
    return convert(T, Rmath.qbinom(q, n, p, false, false))
end
function binominvlogcdf(n::Real, p::Real, lq::Real)
    T = float(Base.promote_typeof(n, p, lq))
    return convert(T, Rmath.qbinom(lq, n, p, true, true))
end
function binominvlogccdf(n::Real, p::Real, lq::Real)
    T = float(Base.promote_typeof(n, p, lq))
    return convert(T, Rmath.qbinom(lq, n, p, false, true))
end
