using StatsFuns
using Base.Test

function check_rmath(fname, params, aname, a, v, rv, isprob, rtol)
    if isprob
        rd = abs(v / rv - 1.0)
        if rd > rtol
            error("$fname deviates too much from Rmath at " *
                "params = $params, $aname = $a:\n" *
                "  v = $v (rv = $rv)\n  |v/rv - 1| = $rd > $rtol.")
        end
    else
        τ = (1.0 + abs(rv)) * rtol
        ad = abs(v - rv)
        if ad > τ
            error("$fname deviates too much from Rmath at " *
                "params = $params, $aname = $a:\n" *
                "  v = $v (rv = $rv)\n  |v - rv| = $ad > $τ.")
        end
    end
end

macro generate_rmath_compfun(basename)
    # function symbols
    compfun = symbol(string("rmathcomp_", basename))

    pdf     = symbol(string(basename, "pdf"))
    logpdf  = symbol(string(basename, "logpdf"))
    cdf     = symbol(string(basename, "cdf"))
    ccdf    = symbol(string(basename, "ccdf"))
    logcdf  = symbol(string(basename, "logcdf"))
    logccdf = symbol(string(basename, "logccdf"))

    invcdf     = symbol(string(basename, "invcdf"))
    invccdf    = symbol(string(basename, "invccdf"))
    invlogcdf  = symbol(string(basename, "invlogcdf"))
    invlogccdf = symbol(string(basename, "invlogccdf"))

    quote
        function $(esc(compfun)){N}(params::NTuple{N,Real}, X::AbstractArray, rtol=1.0e-12)
            for i = 1:length(X)
                x = X[i]
                check_rmath($(string(pdf)), params, "x", x,
                    $(pdf)(params..., x), RMath.$(pdf)(params..., x), true, rtol)
                check_rmath($(string(logpdf)), params, "x", x,
                    $(logpdf)(params..., x), RMath.$(logpdf)(params..., x), false, rtol)

                check_rmath($(string(cdf)), params, "x", x,
                    $(cdf)(params..., x), RMath.$(cdf)(params..., x), true, rtol)
                check_rmath($(string(ccdf)), params, "x", x,
                    $(ccdf)(params..., x), RMath.$(ccdf)(params..., x), true, rtol)
                check_rmath($(string(logcdf)), params, "x", x,
                    $(logcdf)(params..., x), RMath.$(logcdf)(params..., x), false, rtol)
                check_rmath($(string(logccdf)), params, "x", x,
                    $(logccdf)(params..., x), RMath.$(logccdf)(params..., x), true, rtol)

                p = RMath.$(cdf)(params..., x)
                cp = RMath.$(ccdf)(params..., x)
                lp = log(p)
                lcp = log(cp)

                check_rmath($(string(invcdf)), params, "p", p,
                    $(invcdf)(params..., p), RMath.$(invcdf)(params..., p), false, rtol)
                check_rmath($(string(invccdf)), params, "p", cp,
                    $(invccdf)(params..., cp), RMath.$(invccdf)(params..., cp), false, rtol)
                check_rmath($(string(invlogcdf)), params, "lp", lp,
                    $(invlogcdf)(params..., lp), RMath.$(invlogcdf)(params..., lp), false, rtol)
                check_rmath($(string(invlogccdf)), params, "lp", lcp,
                    $(invlogccdf)(params..., lcp), RMath.$(invlogccdf)(params..., lcp), false, rtol)
            end
        end
    end
end

@generate_rmath_compfun norm
rmathcomp_norm((0.0, 1.0), -6.0:0.01:6.0)
rmathcomp_norm((2.0, 1.0), -3.0:0.01:7.0)
rmathcomp_norm((0.0, 0.5), -3.0:0.01:3.0)
