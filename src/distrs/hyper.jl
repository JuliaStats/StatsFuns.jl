# functions related to hyper-geometric distribution

# R implementations
using .RFunctions:
    # hyperpdf,
    # hyperlogpdf,
    # hypercdf,
    # hyperccdf,
    # hyperlogcdf,
    # hyperlogccdf,
    hyperinvcdf,
    hyperinvccdf,
    hyperinvlogcdf,
    hyperinvlogccdf

function hyperpdf(ms, mf, n, x)
    exp(hyperlogpdf(ms, mf, n, x))
end

hyperlogpdf(ms, mf, n, x) = hyperlogpdf(promote(ms,mf,n,x)...)

function hyperlogpdf(ms::T, mf::T, n::T, x::T) where T<:Real
    #loggamma(ms + 1) - loggamma(x + 1) - loggamma(ms - x + 1) +
    #loggamma(mf + 1) - loggamma(n - x + 1) - loggamma(mf - n + x + 1) -
    #loggamma(ms + mf + 1) + loggamma(n + 1) + loggamma(ms + mf - n + 1)
    -log1p(ms) - logbeta(x + 1, ms - x + 1) -
    log1p(mf) - logbeta(n - x + 1, mf - n + x + 1) +
    log1p(ms+mf) + logbeta(n + 1, ms + mf - n + 1)
end


function _hypercdf(ms, mf, n, x, invert)
    N = ms+mf
    mode = fld((n + 1) * (ms + 1), N + 2)
    local result
    if x < mode
        result = hyperpdf(ms, mf, n, x)
        diff = result
        lower_limit = max(0, n+ms-N)
        while x != lower_limit && diff > (invert ? result : 1)*eps()
            diff *= x * (N + x - n - ms) / ((1+n-x)*(1+ms -x))
            result += diff
            x -= 1
        end
    else
        invert = !invert
        x += 1
        result = hyperpdf(ms, mf, n, x)
        diff = result
        upper_limit = min(ms, n)
        while x != upper_limit && diff > (invert ? 1 : result)*eps()
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
end
