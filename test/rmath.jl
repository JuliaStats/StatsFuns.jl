using StatsFuns
using Rmath: Rmath
using Test

include("utils.jl")

function check_rmath(statsfun, rmathfun, params, a, isprob, rtol)
    v = @inferred(statsfun(params..., a))
    rv = @inferred(rmathfun(a, params...))
    @test v isa float(Base.promote_typeof(params..., a))
    return if isprob
        @test v ≈ oftype(v, rv) rtol = rtol nans = true
    else
        @test v ≈ oftype(v, rv) atol = rtol rtol = rtol nans = true
    end
end

function rmathcomp(basename::String, params, X::AbstractArray, rtol = _default_rtol(params, X))
    rbasename = if basename == "nfdist"
        "nf"
    elseif basename == "ntdist"
        "nt"
    elseif basename == "srdist"
        "tukey"
    else
        basename
    end

    if isdefined(Rmath, Symbol(:d, rbasename))
        stats_pdf = getproperty(@__MODULE__, Symbol(basename, :pdf))
        rmath_pdf = let f = getproperty(Rmath, Symbol(:d, rbasename))
            (a, params...) -> f(a, params..., false)
        end
        @testset "pdf with x=$x" for x in X
            check_rmath(
                stats_pdf, rmath_pdf,
                params, x, true, rtol
            )
        end

        stats_logpdf = getproperty(@__MODULE__, Symbol(basename, :logpdf))
        rmath_logpdf = let f = getproperty(Rmath, Symbol(:d, rbasename))
            (a, params...) -> f(a, params..., true)
        end
        @testset "logpdf with x=$x" for x in X
            check_rmath(
                stats_logpdf, rmath_logpdf,
                params, x, false, rtol
            )
        end
    end

    if isdefined(Rmath, Symbol(:p, rbasename))
        stats_cdf = getproperty(@__MODULE__, Symbol(basename, :cdf))
        rmath_cdf = let f = getproperty(Rmath, Symbol(:p, rbasename))
            (a, params...) -> f(a, params..., true, false)
        end
        @testset "cdf with x=$x" for x in X
            check_rmath(stats_cdf, rmath_cdf, params, x, true, rtol)
        end

        stats_ccdf = getproperty(@__MODULE__, Symbol(basename, :ccdf))
        rmath_ccdf = let f = getproperty(Rmath, Symbol(:p, rbasename))
            (a, params...) -> f(a, params..., false, false)
        end
        @testset "ccdf with x=$x" for x in X
            check_rmath(stats_ccdf, rmath_ccdf, params, x, true, rtol)
        end

        stats_logcdf = getproperty(@__MODULE__, Symbol(basename, :logcdf))
        rmath_logcdf = let f = getproperty(Rmath, Symbol(:p, rbasename))
            (a, params...) -> f(a, params..., true, true)
        end
        @testset "logcdf with x=$x" for x in X
            check_rmath(stats_logcdf, rmath_logcdf, params, x, false, rtol)
        end

        stats_logccdf = getproperty(@__MODULE__, Symbol(basename, :logccdf))
        rmath_logccdf = let f = getproperty(Rmath, Symbol(:p, rbasename))
            (a, params...) -> f(a, params..., false, true)
        end
        @testset "logccdf with x=$x" for x in X
            check_rmath(stats_logccdf, rmath_logccdf, params, x, false, rtol)
        end

        #=
        R version of signrank has system variation
        julia> psignrank(18,10,false,true) # windows
        -0.2076393647782445
        julia> psignrank(18,10,false,true) # linux
        -0.20763936477824452
        This slight difference causes test failures for the inverse functions,
        due to a slight shift in the location of the discontinuity.
    
        This also holds true for wilcox.
        =#
        test_inv = (basename != "signrank" && basename != "wilcox") || !Sys.islinux()
        if isdefined(Rmath, Symbol(:q, rbasename)) && test_inv
            stats_invcdf = getproperty(@__MODULE__, Symbol(basename, :invcdf))
            rmath_invcdf = let f = getproperty(Rmath, Symbol(:q, rbasename))
                (a, params...) -> f(a, params..., true, false)
            end
            p = rmath_cdf.(X, params...)
            @testset "invcdf with q=$_p" for _p in p
                check_rmath(stats_invcdf, rmath_invcdf, params, _p, false, rtol)
            end

            stats_invccdf = getproperty(@__MODULE__, Symbol(basename, :invccdf))
            rmath_invccdf = let f = getproperty(Rmath, Symbol(:q, rbasename))
                (a, params...) -> f(a, params..., false, false)
            end
            cp = rmath_ccdf.(X, params...)
            @testset "invccdf with q=$_p" for _p in cp
                check_rmath(stats_invccdf, rmath_invccdf, params, _p, false, rtol)
            end

            stats_invlogcdf = getproperty(@__MODULE__, Symbol(basename, :invlogcdf))
            rmath_invlogcdf = let f = getproperty(Rmath, Symbol(:q, rbasename))
                (a, params...) -> f(a, params..., true, true)
            end
            lp = rmath_logcdf.(X, params...)
            @testset "invlogcdf with log(q)=$_p" for _p in lp
                check_rmath(stats_invlogcdf, rmath_invlogcdf, params, _p, false, rtol)
            end

            stats_invlogccdf = getproperty(@__MODULE__, Symbol(basename, :invlogccdf))
            rmath_invlogccdf = let f = getproperty(Rmath, Symbol(:q, rbasename))
                (a, params...) -> f(a, params..., false, true)
            end
            lcp = rmath_logccdf.(X, params...)
            @testset "invlogccdf with log(q)=$_p" for _p in lcp
                check_rmath(stats_invlogccdf, rmath_invlogccdf, params, _p, false, rtol)
            end
        end
    end

    return nothing
end

function rmathcomp_tests(basename::String, configs)
    return @testset "$basename" begin
        @testset "params: $params" for (params, data) in configs
            rmathcomp(basename, params, data)
        end
    end
end

### Test cases

@testset "RMath" begin
    rmathcomp_tests(
        "beta", [
            ((0.1, 1.0), 0.0:0.01:1.0),
            ((1.0, 1.0), 0.0:0.01:1.0),
            ((2.0, 3.0), 0.0:0.01:1.0),
            ((10.0, 2.0), 0.0:0.01:1.0),
            ((10, 2), 0.0:0.01:1.0),
            ((1.0f0, 1.0f0), 0.0f0:0.01f0:1.0f0),
            ((1.0, 1.0), 0.0f0:0.01f0:1.0f0),
            ((Float16(1), Float16(1)), Float16(0):Float16(0.01):Float16(1)),
            ((1.0f0, 1.0f0), Float16(0):Float16(0.01):Float16(1)),
            ((10, 2), [0, 1]),
            ((10, 2), (0 // 1):(1 // 100):(1 // 1)),
        ]
    )
    # It is not possible to maintain a rtol of 1e-14 for the cdf when there is such a large difference
    # in the magnitude of the parameters. At least not with the current parameters. Furthermore, it
    # seems that Rmath is actually less accurate than we are but since we are comparing against Rmath
    # we have to use rtol=1e-12 although we are probably only off by around 1e-13.
    @testset "param: (1000, 2)" begin
        rmathcomp(
            "beta",
            (1000, 2),
            # We have to drop the 0.48 value since the R quantile function fails while we succeed.
            # It's tested separate below.
            setdiff(collect(0.0:0.01:1.0), 0.48),
            1.0e-12
        )
        # Test p=0.48 separately since R fails. (It's pretty slow, though, caused by the cdf being 9.0797754e-317)
        @test betainvcdf(1000, 2, betacdf(1000, 2, 0.48)) ≈ 0.48
    end
    @testset "$(f)(0, 0, $x) should error" for f in (betacdf, betaccdf, betalogcdf, betalogccdf),
            x in (0.0, 0.5, 1.0)
        @test_throws DomainError f(0.0, 0.0, x)
    end
    # Have to check them separately since Rmath is not completely consistent with our conventions:
    # - betacdf(α, β, x) = P(X ≤ x) where X ~ Beta(α, β)
    # - betaccdf(α, β, x) = P(X > x) where X ~ Beta(α, β)
    # - betainvcdf(α, β, p) = inf { x : p ≤ P(X ≤ x) } where X ~ Beta(α, β)
    # - betainvccdf(α, β, p) = sup { x : p ≤ P(X > x) } where X ~ Beta(α, β)
    @testset "beta: degenerate cases" begin
        # Check degenerate cases Beta(0, β) and Dirac(α, 0) with α, β > 0
        # Beta(0, β) is a Dirac distribution at x=0
        # Beta(α, 0) is a Dirac distribution at x=1
        α = β = 1 // 2

        for x in 0.0f0:0.01f0:1.0f0
            # Check betacdf
            @test @inferred(betacdf(0, β, x)) === 1.0f0
            @test @inferred(betacdf(α, 0, x)) === (x < 1 ? 0.0f0 : 1.0f0)

            # Check betaccdf, betalogcdf, and betalogccdf based on betacdf
            @test @inferred(betaccdf(0, β, x)) === 1 - betacdf(0, β, x)
            @test @inferred(betaccdf(α, 0, x)) === 1 - betacdf(α, 0, x)
            @test @inferred(betalogcdf(0, β, x)) === log(betacdf(0, β, x))
            @test @inferred(betalogcdf(α, 0, x)) === log(betacdf(α, 0, x))
            @test @inferred(betalogccdf(0, β, x)) === log(betaccdf(0, β, x))
            @test @inferred(betalogccdf(α, 0, x)) === log(betaccdf(α, 0, x))
        end

        for p in 0.0f0:0.01f0:1.0f0
            # Check betainvcdf
            @test @inferred(betainvcdf(0, β, p)) === 0.0f0
            @test @inferred(betainvcdf(α, 0, p)) === (p > 0 ? 1.0f0 : 0.0f0)

            # Check betainvccdf
            @test @inferred(betainvccdf(0, β, p)) === (p > 0 ? 0.0f0 : 1.0f0)
            @test @inferred(betainvccdf(α, 0, p)) === 1.0f0
        end
    end

    rmathcomp_tests(
        "binom", [
            ((1, 0.5), 0.0:1.0),
            ((1, 0.7), 0.0:1.0),
            ((8, 0.6), 0.0:8.0),
            ((20, 0.1), 0.0:20.0),
            ((20, 0.9), 0.0:20.0),
            ((20, 0.9), 0:20),
            ((1, 0.5f0), 0.0f0:1.0f0),
            ((1, 0.5), 0.0f0:1.0f0),
            ((1, Float16(0.5)), Float16(0):Float16(1)),
            ((1, 0.5f0), Float16(0):Float16(1)),
            ((10, 1 // 2), (0 // 1):(10 // 1)),
        ]
    )

    rmathcomp_tests(
        "chisq", [
            ((1,), 0.0:0.1:8.0),
            ((4,), 0.0:0.1:8.0),
            ((9,), 0.0:0.1:8.0),
            ((9,), 0:8),
            ((1,), 0.0f0:0.1f0:8.0f0),
            ((1,), Float16(0):Float16(0.1):Float16(8)),
            ((9,), (0 // 1):(8 // 1)),
        ]
    )

    rmathcomp_tests(
        "fdist", [
            ((1, 1), (0.0:0.1:5.0)),
            ((2, 1), (0.0:0.1:5.0)),
            ((5, 2), (0.0:0.1:5.0)),
            ((10, 1), (0.0:0.1:5.0)),
            ((10, 3), (0.0:0.1:5.0)),
            ((10, 3), (0:5)),
            ((1, 1), (0.0f0:0.1f0:5.0f0)),
            ((1, 1), (Float16(0):Float16(0.1):Float16(5))),
            ((10, 3), (0 // 1):(5 // 1)),
        ]
    )

    rmathcomp_tests(
        "gamma", [
            ((1.0, 1.0), (0.0:0.05:12.0)),
            ((0.5, 1.0), (0.0:0.05:12.0)),
            ((3.0, 1.0), (0.0:0.05:12.0)),
            ((9.0, 1.0), (0.0:0.05:12.0)),
            ((2.0, 3.0), (0.0:0.05:12.0)),
            ((2, 3), (0:12)),
            ((1.0f0, 1.0f0), (0.0f0:0.05f0:12.0f0)),
            ((1.0, 1.0), (0.0f0:0.05f0:12.0f0)),
            ((Float16(1), Float16(1)), (Float16(0):Float16(0.05):Float16(12))),
            ((1.0f0, 1.0f0), (Float16(0):Float16(0.05):Float16(12))),
            ((2, 3), ((0 // 1):(12 // 1))),
            ((3.0, 1.0e-310), (0.0:0.05:12.0)),
        ]
    )

    rmathcomp_tests(
        "hyper", [
            ((2, 3, 4), 0.0:4.0),
            ((2, 3, 4), 0:4),
            ((2, 3, 4), 0.0f0:4.0f0),
            ((2, 3, 4), Float16(0):Float16(4)),
            ((2, 3, 4), (0 // 1):(4 // 1)),
        ]
    )

    rmathcomp_tests(
        "nbeta", [
            ((1.0, 1.0, 0.0), 0.01:0.01:0.99),
            ((2.0, 3.0, 0.0), 0.01:0.01:0.99),
            ((1.0, 1.0, 2.0), 0.01:0.01:0.99),
            ((3.0, 4.0, 2.0), 0.01:0.01:0.99),
            ((3, 4, 2), 0.01:0.01:0.99),
            ((1.0f0, 1.0f0, 0.0f0), 0.01f0:0.01f0:0.99f0),
            ((1.0, 1.0, 0.0), 0.01f0:0.01f0:0.99f0),
            ((Float16(1), Float16(1), Float16(0)), Float16(0.01):Float16(0.01):Float16(0.99)),
            ((1.0f0, 1.0f0, 0.0f0), Float16(0.01):Float16(0.01):Float16(0.99)),
            ((3, 4, 2), (1 // 100):(1 // 100):(99 // 100)),
        ]
    )

    rmathcomp_tests(
        "nbinom", [
            ((1, 0.5), 0.0:20.0),
            ((3, 0.5), 0.0:20.0),
            ((3, 0.2), 0.0:20.0),
            ((3, 0.8), 0.0:20.0),
            ((1, 0.5f0), 0.0f0:20.0f0),
            ((1, 0.5), 0.0f0:20.0f0),
            ((1, Float16(0.5)), Float16(0):Float16(20)),
            ((3, 0.8), 0:20),
            ((3, 1 // 2), (0 // 1):(20 // 1)),
        ]
    )

    rmathcomp_tests(
        "nchisq", [
            ((2, 1), 0.0:0.2:8.0),
            ((2, 3), 0.0:0.2:8.0),
            ((4, 1), 0.0:0.2:8.0),
            ((4, 3), 0.0:0.2:8.0),
            ((4, 3), 0:8),
            ((2, 1), 0.0f0:0.2f0:8.0f0),
            ((2, 1), Float16(0):Float16(0.2):Float16(8)),
            ((2, 1), (0 // 1):(1 // 5):(8 // 1)),
        ]
    )

    rmathcomp_tests(
        "nfdist", [
            ((1.0, 1.0, 0.0), 0.1:0.1:10.0),
            ((1.0, 1.0, 2.0), 0.1:0.1:10.0),
            ((2.0, 3.0, 1.0), 0.1:0.1:10.0),
            ((2, 3, 1), 1:10),
            ((1.0f0, 1.0f0, 0.0f0), 0.1f0:0.1f0:10.0f0),
            ((1.0, 1.0, 0.0), 0.1f0:0.1f0:10.0f0),
            ((Float16(1), Float16(1), Float16(0)), Float16(0.1):Float16(0.1):Float16(10)),
            ((1.0f0, 1.0f0, 0.0f0), Float16(0.1):Float16(0.1):Float16(10)),
            ((2, 3, 1), (1 // 1):(10 // 1)),
        ]
    )

    rmathcomp_tests(
        "norm", [
            ((0.0, 1.0), -6.0:0.01:6.0),
            ((2.0, 1.0), -3.0:0.01:7.0),
            ((0.0, 0.5), -3.0:0.01:3.0),
            ((0, 1), -3.0:3.0),
            ((0, 1), -3.0:0.01:3.0),
            ((0, 1), -3:3),
            ((0.0, 0.0), -6.0:0.1:6.0),
            ((0.0f0, 1.0f0), -6.0f0:0.01f0:6.0f0),
            ((0.0, 1.0), -6.0f0:0.01f0:6.0f0),
            ((0, 2), (-6 // 1):(1 // 2):(6 // 1)),
            ((0.0f0, 2.0f0), (-6 // 1):(1 // 2):(6 // 1)),
            # Fail since `SpecialFunctions.erfcx` is not implemented for `Float16`
            #((Float16(0), Float16(1)), -Float16(6):Float16(0.01):Float16(6)),
            #((0f0, 1f0), -Float16(6):Float16(0.01):Float16(6)),
        ]
    )

    rmathcomp_tests(
        "ntdist", [
            ((0, 1), -4.0:0.1:10.0),
            ((0, 4), -4.0:0.1:10.0),
            ((2, 1), -4.0:0.1:10.0),
            ((2, 4), -4.0:0.1:10.0),
            ((2, 4), -4:10),
            ((0, 1), -4.0f0:0.1f0:10.0f0),
            ((0, 1), -Float16(4):Float16(0.1):Float16(10)),
            ((0, 1), (-4 // 1):(1 // 10):(10 // 1)),
        ]
    )

    rmathcomp_tests(
        "pois", [
            ((0.5,), 0.0:20.0),
            ((1.0,), 0.0:20.0),
            ((2.0,), 0.0:20.0),
            ((10.0,), 0.0:20.0),
            ((10,), 0:20),
            ((0.5f0,), 0.0f0:20.0f0),
            ((0.5,), 0.0f0:20.0f0),
            ((Float16(0.5),), Float16(0):Float16(20)),
            ((0.5f0,), Float16(0):Float16(20)),
            ((1 // 2,), (0 // 1):(20 // 1)),
        ]
    )

    rmathcomp_tests(
        "tdist", [
            ((1,), -5.0:0.1:5.0),
            ((2,), -5.0:0.1:5.0),
            ((5,), -5.0:0.1:5.0),
            ((5,), -5:5),
            ((1,), -5.0f0:0.1f0:5.0f0),
            ((1,), -Float16(5):Float16(0.1):Float16(5)),
            ((1,), (-5 // 1):(5 // 1)),
            ((Inf,), -5.0:0.1:5.0),
            ((Inf32,), -5.0f0:0.1f0:5.0f0),
        ]
    )

    rmathcomp_tests(
        "signrank", [
            ((4), -2:12),
            ((4), -2.0:0.25:12.0),
            ((10), -2:57),
        ]
    )

    # The R version does not roundtrip cdf->invcdf, while our version does
    @test_broken RFunctions.signrankinvcdf.(50, RFunctions.signrankcdf.(50, -1:1276)) == [0; 0:1275; 1275]
    @test signrankinvcdf.(50, signrankcdf.(50, -1:1276)) == [0; 0:1275; 1275]
    @test signrankinvccdf.(50, signrankccdf.(50, -1:1276)) == [0; 0:1275; 1275]
    @test signrankinvlogcdf.(50, signranklogcdf.(50, 0:1276)) == [0:1275; 1275]
    @test isnan(signrankinvlogcdf.(50, signranklogcdf(500, -1)))
    @test signrankinvlogccdf.(50, signranklogccdf.(50, -1:1274)) == [0; 0:1274]
    @test isnan(signrankinvlogccdf.(50, signranklogccdf.(50, 1275)))
    @test isnan(signrankinvlogccdf.(50, signranklogccdf.(50, 1276)))

    rmathcomp_tests(
        "srdist", [
            ((1, 2), (0.0:0.2:5.0)),
            ((2, 2), (0.0:0.2:5.0)),
            ((5, 3), (0.0:0.2:5.0)),
            ((10, 2), (0.0:0.2:5.0)),
            ((10, 5), (0.0:0.2:5.0)),
            ((1, 2), (0.0f0:0.2f0:5.0f0)),
            ((1, 2), (Float16(0):Float16(0.2):Float16(5))),
            ((1, 2), ((0 // 1):(1 // 5):(5 // 1))),
        ]
    )

    rmathcomp_tests(
        "wilcox", [
            ((3, 4), -2:13),
            ((3, 4), -2.0:13.0),
            ((4, 3), -2:13),
            ((4, 4), -2:17),
            ((1, 7), -2:8),
            ((1, 7), -2:8),
        ]
    )


    @test wilcoxinvcdf.(10, 10, wilcoxcdf.(10, 10, -1:101)) == [0; 0:100; 100]
    @test wilcoxinvccdf.(10, 10, wilcoxccdf.(10, 10, -1:101)) == [0; 0:100; 100]
    @test wilcoxinvlogcdf.(10, 10, wilcoxlogcdf.(10, 10, 0:101)) == [0:100; 100]
    @test isnan(wilcoxinvlogcdf.(10, 10, wilcoxlogcdf(10, 10, -1)))
    @test wilcoxinvlogccdf.(10, 10, wilcoxlogccdf.(10, 10, -1:99)) == [0; 0:99]
    @test isnan(wilcoxinvlogccdf.(10, 10, wilcoxlogccdf.(10, 10, 100)))
    @test isnan(wilcoxinvlogccdf.(10, 10, wilcoxlogccdf.(10, 10, 101)))

    # Note: Convergence fails in srdist with cdf values below 0.16 with df = 10, k = 5.
    # Reduced df or k allows convergence. This test documents this behavior.
    x = 0.15
    q = srdistcdf(10, 5, 0.15)
    rx = srdistinvcdf(10, 5, q)
    rtol = 100eps(1.0)
    @test_broken x ≈ rx atol = rtol rtol = rtol nans = true

    # Test values outside of the support
    rmathcomp_tests(
        "beta", [
            ((1.0, 1.0), [-10.0, -6.3, 2.1, 23.5]),
            ((1 // 1, 1 // 1), [-10 // 1, -63 // 10, 21 // 10, 47 // 2]),
            ((1, 1), [-10, -6, 2, 24]),
        ]
    )
    rmathcomp_tests(
        "binom", [
            ((5, 0.5), [-8, -2.3, 1.2, 5.4, 11.9]),
            ((5, 1 // 2), [-8, -23 // 10, 6 // 5, 27 // 5, 119 // 10]),
            ((5, 1 // 2), [-8, -2, 6, 12]),
        ]
    )
    rmathcomp_tests(
        "fdist", [
            ((1.0, 1.0), [-10.0, -6.3]),
            ((1 // 1, 1 // 1), [-10 // 1, -63 // 10]),
            ((1, 1), [-10, -6]),
        ]
    )
    rmathcomp_tests(
        "gamma", [
            ((1.0, 1.0), [-10.0, -6.3]),
            ((1 // 1, 1 // 1), [-10 // 1, -63 // 10]),
            ((1, 1), [-10, -6]),
        ]
    )
    rmathcomp_tests(
        "pois", [
            ((0.5,), [-10, -2.5, 1.3, 8.7]),
            ((1 // 2,), [-10, -5 // 2, 13 // 10, 87 // 10]),
            ((1,), [-10, -3]),
        ]
    )
end
