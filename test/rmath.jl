using StatsFuns
using Base.Test
import StatsFuns.Rmath

function check_rmath(fname, statsfun, rmathfun, params, aname, a, isprob, rtol)
    v = statsfun(params..., a)
    rv = rmathfun(params..., a)
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

get_statsfun(fname) = eval(symbol(fname))
get_rmathfun(fname) = eval(parse(string("Rmath.", fname)))

function rmathcomp(basename, params, X::AbstractArray, rtol=1.0e-12)
    pdf     = string(basename, "pdf")
    logpdf  = string(basename, "logpdf")
    cdf     = string(basename, "cdf")
    ccdf    = string(basename, "ccdf")
    logcdf  = string(basename, "logcdf")
    logccdf = string(basename, "logccdf")
    invcdf     = string(basename, "invcdf")
    invccdf    = string(basename, "invccdf")
    invlogcdf  = string(basename, "invlogcdf")
    invlogccdf = string(basename, "invlogccdf")

    stats_pdf     = get_statsfun(pdf)
    stats_logpdf  = get_statsfun(logpdf)
    stats_cdf     = get_statsfun(cdf)
    stats_ccdf    = get_statsfun(ccdf)
    stats_logcdf  = get_statsfun(logcdf)
    stats_logccdf = get_statsfun(logccdf)
    stats_invcdf     = get_statsfun(invcdf)
    stats_invccdf    = get_statsfun(invccdf)
    stats_invlogcdf  = get_statsfun(invlogcdf)
    stats_invlogccdf = get_statsfun(invlogccdf)

    rmath_pdf     = get_rmathfun(pdf)
    rmath_logpdf  = get_rmathfun(logpdf)
    rmath_cdf     = get_rmathfun(cdf)
    rmath_ccdf    = get_rmathfun(ccdf)
    rmath_logcdf  = get_rmathfun(logcdf)
    rmath_logccdf = get_rmathfun(logccdf)
    rmath_invcdf     = get_rmathfun(invcdf)
    rmath_invccdf    = get_rmathfun(invccdf)
    rmath_invlogcdf  = get_rmathfun(invlogcdf)
    rmath_invlogccdf = get_rmathfun(invlogccdf)

    for i = 1:length(X)
        x = X[i]
        check_rmath(pdf, stats_pdf, rmath_pdf,
            params, "x", x, true, rtol)
        check_rmath(logpdf, stats_logpdf, rmath_logpdf,
            params, "x", x, false, rtol)
        check_rmath(cdf, stats_cdf, rmath_cdf,
            params, "x", x, true, rtol)
        check_rmath(ccdf, stats_ccdf, rmath_ccdf,
            params, "x", x, true, rtol)
        check_rmath(logcdf, stats_logcdf, rmath_logcdf,
            params, "x", x, false, rtol)
        check_rmath(logccdf, stats_logccdf, rmath_logccdf,
            params, "x", x, false, rtol)

        p = rmath_cdf(params..., x)
        cp = rmath_ccdf(params..., x)
        lp = rmath_logcdf(params..., x)
        lcp = rmath_logccdf(params..., x)

        check_rmath(invcdf, stats_invcdf, rmath_invcdf,
            params, "q", p, false, rtol)
        check_rmath(invccdf, stats_invccdf, rmath_invccdf,
            params, "q", cp, false, rtol)
        check_rmath(invlogcdf, stats_invlogcdf, rmath_invlogcdf,
            params, "lq", lp, false, rtol)
        check_rmath(invlogccdf, stats_invlogccdf, rmath_invlogccdf,
            params, "lq", lcp, false, rtol)
    end
end

### Test cases

println("\ttesting beta ...")
rmathcomp("beta", (1.0, 1.0), 0.01:0.01:0.99)
rmathcomp("beta", (2.0, 3.0), 0.01:0.01:0.99)
rmathcomp("beta", (10.0, 2.0), 0.01:0.01:0.99)

println("\ttesting binom ...")
rmathcomp("binom", (1, 0.5), 0:1)
rmathcomp("binom", (1, 0.7), 0:1)
rmathcomp("binom", (8, 0.6), 0:8)
rmathcomp("binom", (20, 0.1), 0:20)
rmathcomp("binom", (20, 0.9), 0:20)

println("\ttesting norm ...")
rmathcomp("norm", (0.0, 1.0), -6.0:0.01:6.0)
rmathcomp("norm", (2.0, 1.0), -3.0:0.01:7.0)
rmathcomp("norm", (0.0, 0.5), -3.0:0.01:3.0)
