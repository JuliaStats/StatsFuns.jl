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
    rand       = string(basename, "rand")

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

    # tackle rand specially
    has_rand = true
    if basename == "nbeta" || basename == "nfdist" || basename == "ntdist"
        has_rand = false
    end
    rmath_rand = has_rand ? get_rmathfun(rand) : nothing

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

        # make sure that rand works
        if has_rand
            rmath_rand(params...)
        end
    end
end

function rmathcomp_tests(basename, configs)
    println("\ttesting $basename ...")
    for (params, data) in configs
        rmathcomp(basename, params, data)
    end
end

### Test cases

rmathcomp_tests("beta", [
    ((1.0, 1.0), 0.01:0.01:0.99),
    ((2.0, 3.0), 0.01:0.01:0.99),
    ((10.0, 2.0), 0.01:0.01:0.99),
])

rmathcomp_tests("binom", [
    ((1, 0.5), 0:1),
    ((1, 0.7), 0:1),
    ((8, 0.6), 0:8),
    ((20, 0.1), 0:20),
    ((20, 0.9), 0:20),
])

rmathcomp_tests("chisq", [
    ((1,), 0.0:0.1:8.0),
    ((4,), 0.0:0.1:8.0),
    ((9,), 0.0:0.1:8.0),
])

rmathcomp_tests("fdist", [
    ((1, 1), (0.0:0.1:5.0)),
    ((2, 1), (0.0:0.1:5.0)),
    ((5, 2), (0.0:0.1:5.0)),
    ((10, 1), (0.0:0.1:5.0)),
    ((10, 3), (0.0:0.1:5.0)),
])

rmathcomp_tests("gamma", [
    ((1.0, 1.0), (0.05:0.05:12.0)),
    ((0.5, 1.0), (0.05:0.05:12.0)),
    ((3.0, 1.0), (0.05:0.05:12.0)),
    ((9.0, 1.0), (0.05:0.05:12.0)),
    ((2.0, 3.0), (0.05:0.05:12.0)),
])

rmathcomp_tests("hyper", [
    ((2, 3, 4), 0:4)
])

rmathcomp_tests("nbeta", [
    ((1.0, 1.0, 0.0), 0.01:0.01:0.99),
    ((2.0, 3.0, 0.0), 0.01:0.01:0.99),
    ((1.0, 1.0, 2.0), 0.01:0.01:0.99),
    ((3.0, 4.0, 2.0), 0.01:0.01:0.99),
])

rmathcomp_tests("nbinom", [
    ((1, 0.5), 0:20),
    ((3, 0.5), 0:20),
    ((3, 0.2), 0:20),
    ((3, 0.8), 0:20),
])

rmathcomp_tests("nchisq", [
    ((2, 1), 0.0:0.2:8.0),
    ((2, 3), 0.0:0.2:8.0),
    ((4, 1), 0.0:0.2:8.0),
    ((4, 3), 0.0:0.2:8.0),
])

rmathcomp_tests("nfdist", [
    ((1.0, 1.0, 0.0), 0.1:0.1:10.0),
    ((1.0, 1.0, 2.0), 0.1:0.1:10.0),
    ((2.0, 3.0, 1.0), 0.1:0.1:10.0),
])

rmathcomp_tests("norm", [
    ((0.0, 1.0), -6.0:0.01:6.0),
    ((2.0, 1.0), -3.0:0.01:7.0),
    ((0.0, 0.5), -3.0:0.01:3.0),
])

rmathcomp_tests("ntdist", [
    ((0, 1), -4.0:0.1:10.0),
    ((0, 4), -4.0:0.1:10.0),
    ((2, 1), -4.0:0.1:10.0),
    ((2, 4), -4.0:0.1:10.0),
])

rmathcomp_tests("pois", [
    ((0.5,), 0:20),
    ((1.0,), 0:20),
    ((2.0,), 0:20),
    ((10.0,), 0:20),
])

rmathcomp_tests("tdist", [
    ((1,), -5.0:0.1:5.0),
    ((2,), -5.0:0.1:5.0),
    ((5,), -5.0:0.1:5.0),
])
