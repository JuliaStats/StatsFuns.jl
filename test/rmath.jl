using StatsFuns, Test
import StatsFuns.RFunctions

function check_rmath(fname, statsfun, rmathfun, params, aname, a, isprob, rtol)
    v = @inferred(statsfun(params..., a))
    rv = @inferred(rmathfun(params..., a))
    @test rv isa float(Base.promote_typeof(params..., a))
    if isprob
        @test v ≈ rv rtol=rtol nans=true
    else
        @test v ≈ rv atol=rtol rtol=rtol nans=true
    end
end

get_statsfun(fname) = eval(Symbol(fname))
get_rmathfun(fname) = eval(Meta.parse(string("RFunctions.", fname)))

function rmathcomp(basename, params, X::AbstractArray, rtol=100eps(float(one(eltype(X)))))
    # tackle pdf specially
    has_pdf = true
    if basename == "srdist"
        has_pdf = false
    end

    if has_pdf
        pdf     = string(basename, "pdf")
        logpdf  = string(basename, "logpdf")
    end
    cdf     = string(basename, "cdf")
    ccdf    = string(basename, "ccdf")
    logcdf  = string(basename, "logcdf")
    logccdf = string(basename, "logccdf")
    invcdf     = string(basename, "invcdf")
    invccdf    = string(basename, "invccdf")
    invlogcdf  = string(basename, "invlogcdf")
    invlogccdf = string(basename, "invlogccdf")
    rand       = string(basename, "rand")

    if has_pdf
        stats_pdf     = get_statsfun(pdf)
        stats_logpdf  = get_statsfun(logpdf)
    end
    stats_cdf     = get_statsfun(cdf)
    stats_ccdf    = get_statsfun(ccdf)
    stats_logcdf  = get_statsfun(logcdf)
    stats_logccdf = get_statsfun(logccdf)
    stats_invcdf     = get_statsfun(invcdf)
    stats_invccdf    = get_statsfun(invccdf)
    stats_invlogcdf  = get_statsfun(invlogcdf)
    stats_invlogccdf = get_statsfun(invlogccdf)

    if has_pdf
        rmath_pdf     = get_rmathfun(pdf)
        rmath_logpdf  = get_rmathfun(logpdf)
    end
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
    if basename == "nbeta" || basename == "nfdist" || basename == "ntdist" || basename == "srdist"
        has_rand = false
    end
    rmath_rand = has_rand ? get_rmathfun(rand) : nothing

    for x in X
        if has_pdf
            check_rmath(pdf, stats_pdf, rmath_pdf,
                params, "x", x, true, rtol)
            check_rmath(logpdf, stats_logpdf, rmath_logpdf,
                params, "x", x, false, rtol)
        end
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

@testset "RMath" begin
    rmathcomp_tests("beta", [
        ((1.0, 1.0), 0.01:0.01:0.99),
        ((2.0, 3.0), 0.01:0.01:0.99),
        ((10.0, 2.0), 0.01:0.01:0.99),
        ((10, 2), 0.01:0.01:0.99),
        ((1f0, 1f0), 0.01f0:0.01f0:0.99f0),
        ((1.0, 1.0), 0.01f0:0.01f0:0.99f0),
        ((Float16(1), Float16(1)), Float16(0.01):Float16(0.01):Float16(0.99)),
        ((1f0, 1f0), Float16(0.01):Float16(0.01):Float16(0.99)),
    ])

    rmathcomp_tests("binom", [
        ((1, 0.5), 0.0:1.0),
        ((1, 0.7), 0.0:1.0),
        ((8, 0.6), 0.0:8.0),
        ((20, 0.1), 0.0:20.0),
        ((20, 0.9), 0.0:20.0),
        ((20, 0.9), 0:20),
        ((1, 0.5f0), 0f0:1f0),
        ((1, 0.5), 0f0:1f0),
        ((1, Float16(0.5)), Float16(0):Float16(1)),
        ((1, 0.5f0), Float16(0):Float16(1)),
    ])

    rmathcomp_tests("chisq", [
        ((1,), 0.0:0.1:8.0),
        ((4,), 0.0:0.1:8.0),
        ((9,), 0.0:0.1:8.0),
        ((9,), 0:8),
        ((1,), 0f0:0.1f0:8f0),
        ((1,), Float16(0):Float16(0.1):Float16(8)),
    ])

    rmathcomp_tests("fdist", [
        ((1, 1), (0.0:0.1:5.0)),
        ((2, 1), (0.0:0.1:5.0)),
        ((5, 2), (0.0:0.1:5.0)),
        ((10, 1), (0.0:0.1:5.0)),
        ((10, 3), (0.0:0.1:5.0)),
        ((10, 3), (0:5)),
        ((1, 1), (0f0:0.1f0:5f0)),
        ((1, 1), (Float16(0):Float16(0.1):Float16(5))),
    ])

    rmathcomp_tests("gamma", [
        ((1.0, 1.0), (0.05:0.05:12.0)),
        ((0.5, 1.0), (0.05:0.05:12.0)),
        ((3.0, 1.0), (0.05:0.05:12.0)),
        ((9.0, 1.0), (0.05:0.05:12.0)),
        ((2.0, 3.0), (0.05:0.05:12.0)),
        ((2, 3), (1:12)),
        ((1f0, 1f0), (0.05f0:0.05f0:12f0)),
        ((1.0, 1.0), (0.05f0:0.05f0:12f0)),
        ((Float16(1), Float16(1)), (Float16(0.05):Float16(0.05):Float16(12))),
        ((1f0, 1f0), (Float16(0.05):Float16(0.05):Float16(12))),
    ])

    rmathcomp_tests("hyper", [
        ((2, 3, 4), 0.0:4.0),
        ((2, 3, 4), 0:4),
        ((2, 3, 4), 0f0:4f0),
        ((2, 3, 4), Float16(0):Float16(4)),
    ])

    rmathcomp_tests("nbeta", [
        ((1.0, 1.0, 0.0), 0.01:0.01:0.99),
        ((2.0, 3.0, 0.0), 0.01:0.01:0.99),
        ((1.0, 1.0, 2.0), 0.01:0.01:0.99),
        ((3.0, 4.0, 2.0), 0.01:0.01:0.99),
        ((3, 4, 2), 0.01:0.01:0.99),
        ((1f0, 1f0, 0f0), 0.01f0:0.01f0:0.99f0),
        ((1.0, 1.0, 0.0), 0.01f0:0.01f0:0.99f0),
        ((Float16(1), Float16(1), Float16(0)), Float16(0.01):Float16(0.01):Float16(0.99)),
        ((1f0, 1f0, 0f0), Float16(0.01):Float16(0.01):Float16(0.99)),
    ])

    rmathcomp_tests("nbinom", [
        ((1, 0.5), 0.0:20.0),
        ((3, 0.5), 0.0:20.0),
        ((3, 0.2), 0.0:20.0),
        ((3, 0.8), 0.0:20.0),
        ((1, 0.5f0), 0f0:20f0),
        ((1, 0.5), 0f0:20f0),
        ((1, Float16(0.5)), Float16(0):Float16(20)),
        ((3, 0.8), 0:20),
    ])

    rmathcomp_tests("nchisq", [
        ((2, 1), 0.0:0.2:8.0),
        ((2, 3), 0.0:0.2:8.0),
        ((4, 1), 0.0:0.2:8.0),
        ((4, 3), 0.0:0.2:8.0),
        ((4, 3), 0:8),
        ((2, 1), 0f0:0.2f0:8f0),
        ((2, 1), Float16(0):Float16(0.2):Float16(8)),
    ])

    rmathcomp_tests("nfdist", [
        ((1.0, 1.0, 0.0), 0.1:0.1:10.0),
        ((1.0, 1.0, 2.0), 0.1:0.1:10.0),
        ((2.0, 3.0, 1.0), 0.1:0.1:10.0),
        ((2, 3, 1), 1:10),
        ((1f0, 1f0, 0f0), 0.1f0:0.1f0:10f0),
        ((1.0, 1.0, 0.0), 0.1f0:0.1f0:10f0),
        ((Float16(1), Float16(1), Float16(0)), Float16(0.1):Float16(0.1):Float16(10)),
        ((1f0, 1f0, 0f0), Float16(0.1):Float16(0.1):Float16(10)),
    ])

    rmathcomp_tests("norm", [
        ((0.0, 1.0), -6.0:0.01:6.0),
        ((2.0, 1.0), -3.0:0.01:7.0),
        ((0.0, 0.5), -3.0:0.01:3.0),
        ((0, 1), -3.0:3.0),
        ((0, 1), -3.0:0.01:3.0),
        ((0, 1), -3:3),
        ((0.0, 0.0), -6.0:0.1:6.0),
        ((0f0, 1f0), -6f0:0.01f0:6f0),
        ((0.0, 1.0), -6f0:0.01f0:6f0),
        # Fail since `SpecialFunctions.erfcx` is not implemented for `Float16`
        #((Float16(0), Float16(1)), -Float16(6):Float16(0.01):Float16(6)),
        #((0f0, 1f0), -Float16(6):Float16(0.01):Float16(6)),
    ])

    rmathcomp_tests("ntdist", [
        ((0, 1), -4.0:0.1:10.0),
        ((0, 4), -4.0:0.1:10.0),
        ((2, 1), -4.0:0.1:10.0),
        ((2, 4), -4.0:0.1:10.0),
        ((2, 4), -4:10),
        ((0, 1), -4f0:0.1f0:10f0),
        ((0, 1), -Float16(4):Float16(0.1):Float16(10)),
    ])

    rmathcomp_tests("pois", [
        ((0.5,), 0.0:20.0),
        ((1.0,), 0.0:20.0),
        ((2.0,), 0.0:20.0),
        ((10.0,), 0.0:20.0),
        ((10,), 0:20),
        ((0.5f0,), 0f0:20f0),
        ((0.5,), 0f0:20f0),
        ((Float16(0.5),), Float16(0):Float16(20)),
        ((0.5f0,), Float16(0):Float16(20)),
    ])

    rmathcomp_tests("tdist", [
        ((1,), -5.0:0.1:5.0),
        ((2,), -5.0:0.1:5.0),
        ((5,), -5.0:0.1:5.0),
        ((5,), -5:5),
        ((1,), -5f0:0.1f0:5f0),
        ((1,), -Float16(5):Float16(0.1):Float16(5)),
    ])

    rmathcomp_tests("srdist", [
        ((1,2), (0.0:0.2:5.0)),
        ((2,2), (0.0:0.2:5.0)),
        ((5,3), (0.0:0.2:5.0)),
        ((10,2), (0.0:0.2:5.0)),
        ((10,5), (0.0:0.2:5.0)),
        ((1,2), (0f0:0.2f0:5f0)),
        ((1,2), (Float16(0):Float16(0.2):Float16(5))),
    ])

    # Note: Convergence fails in srdist with cdf values below 0.16 with df = 10, k = 5.
    # Reduced df or k allows convergence. This test documents this behavior.
    x = 0.15
    q = srdistcdf(10, 5, 0.15)
    rx = srdistinvcdf(10, 5, q)
    rtol = 100eps(1.0)
    @test_broken x ≈ rx atol=rtol rtol=rtol nans=true
end
