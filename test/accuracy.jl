@testmodule FDistRef begin
    using SpecialFunctions: logbeta

    # References for implementations that have to cope with extreme parameters. These are the
    # textbook formulas, which cancel badly in `Float64` - that is precisely what is being
    # tested, and `BigFloat` has enough precision (256 bits by default) to absorb it.
    function logpdf(ν1::Real, ν2::Real, x::Real)
        _ν1, _ν2, _x = BigFloat(ν1), BigFloat(ν2), BigFloat(x)
        a, b = _ν1 / 2, _ν2 / 2
        return a * log(_ν1 / _ν2) + (a - 1) * log(_x) - (a + b) * log1p(_ν1 * _x / _ν2) - logbeta(a, b)
    end

    # For `x` so small that the beta variate `u = ν1 * x / (ν1 * x + ν2)` is negligible the cdf
    # is `u^a / (a * beta(a, b)) * (1 + O(u))`, which needs no incomplete beta function.
    function logcdf_smallx(ν1::Real, ν2::Real, x::Real)
        _ν1, _ν2, _x = BigFloat(ν1), BigFloat(ν2), BigFloat(x)
        a, b = _ν1 / 2, _ν2 / 2
        u = _ν1 * _x / (_ν1 * _x + _ν2)
        return a * log(u) - log(a) - logbeta(a, b)
    end
end

@testitem "Accuracy" setup = [FDistRef] begin
    using StatsFuns
    using Test

    # the density has to stay accurate when either degree of freedom is huge: for large
    # `ν2` the term `ν1 * x / ν2` is tiny and must not be absorbed by an explicit `1 +`,
    # and for large `ν1` no two terms growing like `ν1 * log(ν1)` may cancel
    @testset "fdistlogpdf with ν$i = $ν" for (i, ν) in
        Iterators.product(1:2, (1.0e9, 1.0e12, 1.0e15, 1.0e18, 1.0e20))

        ν1, ν2 = i == 1 ? (ν, 5.0) : (1.0, ν)
        @testset "x = $x" for x in (1.0e-300, 1.0e-8, 0.5, 2.0, 6.6, 1.0e8)
            ref = FDistRef.logpdf(ν1, ν2, x)
            @test fdistlogpdf(ν1, ν2, x) ≈ ref rtol = 1.0e-12
            # `exp` turns the absolute error of the log-density into a relative one, and the
            # reference has to be rounded first since `BigFloat` neither underflows nor
            # overflows where `Float64` does
            @test fdistpdf(ν1, ν2, x) ≈ Float64(exp(ref)) rtol = 1.0e-12 * max(1, abs(Float64(ref)))
        end
    end
end

# `ν2 / (ν1 * x)` overflows for `x < ν2 / (ν1 * floatmax(T))` - an ordinary argument in
# `Float16` - and used to take the density and the cdf with it
@testitem "fdist with small x" setup = [FDistRef] begin
    using StatsFuns
    using Test

    @testset "T = $T" for (T, x) in ((Float64, 1.0e-300), (Float32, 1.0f-38), (Float16, Float16(1.0e-4)))
        ν1, ν2 = one(T), T(7)

        @test fdistlogpdf(ν1, ν2, x) ≈ FDistRef.logpdf(ν1, ν2, x) rtol = eps(T)^(3 // 4)

        ref = FDistRef.logcdf_smallx(ν1, ν2, x)
        @test fdistlogcdf(ν1, ν2, x) ≈ ref rtol = eps(T)^(2 // 3)
        @test fdistcdf(ν1, ν2, x) ≈ T(exp(ref)) rtol = eps(T)^(2 // 3) * max(1, abs(Float64(ref)))
    end

    # the same cliff sits inside the large-`ν2` regime: `x = 1.0e-290` is still fine
    @test fdistlogpdf(1.0, 1.0e18, 1.0e-291) ≈ FDistRef.logpdf(1.0, 1.0e18, 1.0e-291) rtol = 1.0e-12
end

@testitem "fdistlogpdf boundary values" begin
    using StatsFuns
    using Test

    # outside the support, where the beta variate `ν1 * x / (ν1 * x + ν2)` is not defined
    @testset "x < 0" begin
        @test fdistlogpdf(3, 7, -1.0) == -Inf
        @test fdistlogpdf(3, 7, -Inf) == -Inf
        @test fdistpdf(3, 7, -1.0) == 0

        # `oftype` keeps the branch from widening to `Float64`
        @test fdistlogpdf(3, 7, -1.0f0) === -Inf32
        @test fdistlogpdf(3, 7, Float16(-1)) === Float16(-Inf)
    end

    # `xlog1py(ν2 / 2, ν1 * x / ν2)` is `Inf` here, whereas the previous formulation
    # subtracted two infinities and returned `NaN`
    @testset "x = Inf" begin
        @test fdistlogpdf(3, 7, Inf) == -Inf
        @test fdistpdf(3, 7, Inf) == 0
    end

    # at zero the density behaves like `x^(ν1 / 2 - 1)`; `-0.0` is not `< 0` and takes the
    # same branch as `0.0`
    @testset "x = 0" begin
        @test fdistlogpdf(1, 1, 0.0) == Inf
        @test fdistlogpdf(2, 1, 0.0) == 0
        @test fdistlogpdf(10, 3, 0.0) == -Inf
        @test fdistlogpdf(10, 3, -0.0) == -Inf
    end

    # `NaN` propagates through the `xlogy` of the branch it is sorted into
    @testset "x = NaN" begin
        @test isnan(fdistlogpdf(3, 7, NaN))
        @test isnan(fdistpdf(3, 7, NaN))
    end
end
