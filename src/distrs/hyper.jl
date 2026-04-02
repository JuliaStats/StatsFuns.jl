# functions related to hyper-geometric distribution

hyperpdf(ms::Real, mf::Real, n::Real, x::Real) = hyperpdf(promote(ms,mf,n,x)...)
    
function hyperpdf(ms::T, mf::T, n::T, x::T) where T<:Real
    M, N = ms, ms+mf
    fxnprev = fMNprev = one(T)
    fxn = 1 - 2x/n
    fMN = 1 - 2M/N
    eiprev = ones(T, Int(n+1))
    k = 0:n
    ei = 1 .- k .*(2/n)
    square(x) = x^2
    res = fma(fxn, fMN/sum(square, ei), 1/(n+1))
    ei = collect(ei)
    for i in 1:n-1
        fxnnext = ((2i + 1) * (n - 2x) * fxn - i * (n + i + 1) * fxnprev) / ((i + 1) * (n - i))
        fxn, fxnprev = fxnnext, fxn
        fMNnext = ((2i + 1) * (N - 2M) * fMN - i * (N + i + 1) * fMNprev) / ((i + 1) * (N - i))
        fMN, fMNprev = fMNnext, fMN
        @. eiprev = ((2i + 1) * (n - 2k) * ei - i * (n + i + 1) * eiprev) / ((i + 1) * (n - i))
        ei, eiprev = eiprev, ei
        correction =  fxn*fMN/sum(square, ei)
        abs(correction) < res*eps(T) && break
        res += correction
    end
    res
end

hyperlogpdf(ms::Real, mf::Real, n::Real, x::Real) = hyperlogpdf(promote(ms,mf,n,x)...)
function hyperlogpdf(ms::T, mf::T, n::T, x::T) where T<:Real
    return log1p(ms+mf) - log1p(ms) - log1p(mf) + 
    (logbeta(n + 1, ms + mf - n + 1) - logbeta(x + 1, ms - x + 1) - 
    logbeta(n - x + 1, mf - n + x + 1))
end

_hypercdf(ms::Real, mf::Real, n::Real, x::Real, invert::Bool) = _hypercdf(promote(ms, mf, n, x)..., invert::Bool)
function _hypercdf(ms::T, mf::T, n::T, x::T, invert::Bool) where T<:Real
    N = ms+mf
    mode = fld((n + 1) * (ms + 1), N + 2)
    local result
    if x < mode
        result = diff = hyperpdf(ms, mf, n, x)
        lower_limit = max(0, n+ms-N)
        while x != lower_limit && diff > (invert ? result : 1)*eps(T)
            diff *= x * (N + x - n - ms) / ((1+n-x)*(1+ms -x))
            result += diff
            x -= 1
        end
    else
        invert = !invert
        x += 1
        result = diff = hyperpdf(ms, mf, n, x)
        upper_limit = min(ms, n)
        while x != upper_limit && diff > (invert ? 1 : result)*eps(T)
            diff *= (n-x) * (ms-x) / ((1+x) * (N+x+1-n-ms))
            result += diff
            x += 1
        end
    end
    result, invert
end

function hypercdf(ms, mf, n, x)
    result, invert = _hypercdf(ms, mf, n, x, false)
    invert ? 1-result : result
end

function hyperccdf(ms, mf, n, x)
    result, invert = _hypercdf(ms, mf, n, x, true)
    invert ? 1-result : result
end

function hyperlogcdf(ms, mf, n, x)
    result, invert = _hypercdf(ms, mf, n, x, false)
    invert ? log1p(-result) : log(result)
end
function hyperlogccdf(ms, mf, n, x)
    result, invert = _hypercdf(ms, mf, n, x, true)
    invert ? log1p(-result) : log(result)

# Rmath implementations
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
