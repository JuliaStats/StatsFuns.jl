@testmodule Utils begin
    using StatsFuns
    using Rmath: Rmath
    using Test

    # default relative tolerance for comparisons with Rmath
    function _default_rtol(params, X::AbstractArray)
        # has to take into account `params` as well since otherwise e.g. `X::Array{<:Rational}`
        # always uses a tolerance based on `eps(one(Float64))` even when parameters are of type
        # Float32
        return _default_rtol(float(promote_type(Base.promote_typeof(params...), eltype(X))))
    end

    # We use less sharp tolerances for Float16 and Float32 since the proportion of significant
    # digits that are equal when evaluating Rmath and StatsFuns is smaller
    # eps^(7//8) means requiring equality of about 7/8 of the significant digits, etc.
    _default_rtol(::Type{Float64}) = eps(Float64)^(7 // 8) # ~2.0e-14
    _default_rtol(::Type{Float32}) = eps(Float32)^(3 // 4) # ~6.4e-6
    _default_rtol(::Type{Float16}) = eps(Float16)^(2 // 3) # ~9.9e-3

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

        # Rmath.ptukey / Rmath.qtukey take an extra `nranges` positional argument
        # (defaulting to 1.0) between the distribution parameters and the
        # `lower_tail`/`log_p` flags. Inject it explicitly for tukey so the test
        # lambdas line up with the Rmath C signature.
        extra_rmath_args = rbasename == "tukey" ? (1,) : ()

        if isdefined(Rmath, Symbol(:d, rbasename))
            stats_pdf = getproperty(StatsFuns, Symbol(basename, :pdf))
            rmath_pdf = let f = getproperty(Rmath, Symbol(:d, rbasename))
                (a, params...) -> f(a, params..., false)
            end
            @testset "pdf with x=$x" for x in X
                check_rmath(
                    stats_pdf, rmath_pdf,
                    params, x, true, rtol
                )
            end

            stats_logpdf = getproperty(StatsFuns, Symbol(basename, :logpdf))
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
            stats_cdf = getproperty(StatsFuns, Symbol(basename, :cdf))
            rmath_cdf = let f = getproperty(Rmath, Symbol(:p, rbasename)), extra = extra_rmath_args
                (a, params...) -> f(a, params..., extra..., true, false)
            end
            @testset "cdf with x=$x" for x in X
                check_rmath(stats_cdf, rmath_cdf, params, x, true, rtol)
            end

            stats_ccdf = getproperty(StatsFuns, Symbol(basename, :ccdf))
            rmath_ccdf = let f = getproperty(Rmath, Symbol(:p, rbasename)), extra = extra_rmath_args
                (a, params...) -> f(a, params..., extra..., false, false)
            end
            @testset "ccdf with x=$x" for x in X
                check_rmath(stats_ccdf, rmath_ccdf, params, x, true, rtol)
            end

            stats_logcdf = getproperty(StatsFuns, Symbol(basename, :logcdf))
            rmath_logcdf = let f = getproperty(Rmath, Symbol(:p, rbasename)), extra = extra_rmath_args
                (a, params...) -> f(a, params..., extra..., true, true)
            end
            @testset "logcdf with x=$x" for x in X
                check_rmath(stats_logcdf, rmath_logcdf, params, x, false, rtol)
            end

            stats_logccdf = getproperty(StatsFuns, Symbol(basename, :logccdf))
            rmath_logccdf = let f = getproperty(Rmath, Symbol(:p, rbasename)), extra = extra_rmath_args
                (a, params...) -> f(a, params..., extra..., false, true)
            end
            @testset "logccdf with x=$x" for x in X
                check_rmath(stats_logccdf, rmath_logccdf, params, x, false, rtol)
            end

            #=
            signrank and wilcox are implemented natively rather than by delegating to Rmath, so
            their inverse functions are not expected to reproduce the Rmath quantile functions
            bit for bit. The inverses locate a discontinuity by searching for the first argument
            at which the cdf reaches `q`, and the `q` tested here are themselves cdf values, so
            the comparison is decided by the last bit of the cdf and is not stable across
            platforms. R's own signrank already varies this way:
            julia> psignrank(18,10,false,true) # windows
            -0.2076393647782445
            julia> psignrank(18,10,false,true) # linux
            -0.20763936477824452
            As for the other natively implemented discrete distributions, the inverse functions
            are instead tested by round tripping through our own cdf further down.
            =#
            test_inv = basename != "signrank" && basename != "wilcox"
            if isdefined(Rmath, Symbol(:q, rbasename)) && test_inv
                stats_invcdf = getproperty(StatsFuns, Symbol(basename, :invcdf))
                rmath_invcdf = let f = getproperty(Rmath, Symbol(:q, rbasename)), extra = extra_rmath_args
                    (a, params...) -> f(a, params..., extra..., true, false)
                end
                p = rmath_cdf.(X, params...)
                @testset "invcdf with q=$_p" for _p in p
                    check_rmath(stats_invcdf, rmath_invcdf, params, _p, false, rtol)
                end

                stats_invccdf = getproperty(StatsFuns, Symbol(basename, :invccdf))
                rmath_invccdf = let f = getproperty(Rmath, Symbol(:q, rbasename)), extra = extra_rmath_args
                    (a, params...) -> f(a, params..., extra..., false, false)
                end
                cp = rmath_ccdf.(X, params...)
                @testset "invccdf with q=$_p" for _p in cp
                    check_rmath(stats_invccdf, rmath_invccdf, params, _p, false, rtol)
                end

                stats_invlogcdf = getproperty(StatsFuns, Symbol(basename, :invlogcdf))
                rmath_invlogcdf = let f = getproperty(Rmath, Symbol(:q, rbasename)), extra = extra_rmath_args
                    (a, params...) -> f(a, params..., extra..., true, true)
                end
                lp = rmath_logcdf.(X, params...)
                @testset "invlogcdf with log(q)=$_p" for _p in lp
                    check_rmath(stats_invlogcdf, rmath_invlogcdf, params, _p, false, rtol)
                end

                stats_invlogccdf = getproperty(StatsFuns, Symbol(basename, :invlogccdf))
                rmath_invlogccdf = let f = getproperty(Rmath, Symbol(:q, rbasename)), extra = extra_rmath_args
                    (a, params...) -> f(a, params..., extra..., false, true)
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
end

@testmodule FDistRef begin
    using SpecialFunctions: logbeta

    # textbook formulas for `fdist`, evaluated in `BigFloat`. They cancel badly in `Float64` -
    # that is precisely what is being tested - and for degrees of freedom of `floatmax` scale the
    # cancellation runs to 300 digits, well past the 77 that the default 256 bits give. At 256 bits
    # `logpdf(1e200, 1e300, 0.5)` is off by a relative 3.6e3 and would validate anything
    function logpdf(ν1::Real, ν2::Real, x::Real)
        return setprecision(BigFloat, 1024) do
            _ν1, _ν2, _x = BigFloat(ν1), BigFloat(ν2), BigFloat(x)
            a, b = _ν1 / 2, _ν2 / 2
            a * log(_ν1 / _ν2) + (a - 1) * log(_x) - (a + b) * log1p(_ν1 * _x / _ν2) - logbeta(a, b)
        end
    end

    function dlogpdf_dx(ν1::Real, ν2::Real, x::Real)
        return setprecision(BigFloat, 1024) do
            _ν1, _ν2, _x = BigFloat(ν1), BigFloat(ν2), BigFloat(x)
            (_ν1 / 2 - 1) / _x - ((_ν1 + _ν2) / 2) * _ν1 / (_ν1 * _x + _ν2)
        end
    end

    # for `x` so small that the beta variate `u = ν1 * x / (ν1 * x + ν2)` is negligible the cdf
    # is `u^a / (a * beta(a, b)) * (1 + O(u))`, which needs no incomplete beta function
    function logcdf_smallx(ν1::Real, ν2::Real, x::Real)
        return setprecision(BigFloat, 1024) do
            _ν1, _ν2, _x = BigFloat(ν1), BigFloat(ν2), BigFloat(x)
            a, b = _ν1 / 2, _ν2 / 2
            u = _ν1 * _x / (_ν1 * _x + _ν2)
            a * log(u) - log(a) - logbeta(a, b)
        end
    end
end
