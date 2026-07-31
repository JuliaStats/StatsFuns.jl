@testitem "Accuracy" begin
    using StatsFuns
    using SpecialFunctions: logbeta
    using Test

    # References for implementations that have to cope with extreme parameters. These are the
    # textbook formulas, which cancel badly in `Float64` - that is precisely what is being
    # tested, and `BigFloat` has enough precision (256 bits by default) to absorb it.
    function fdistlogpdf_ref(ν1::Real, ν2::Real, x::Real)
        _ν1, _ν2, _x = BigFloat(ν1), BigFloat(ν2), BigFloat(x)
        a, b = _ν1 / 2, _ν2 / 2
        return a * log(_ν1 / _ν2) + (a - 1) * log(_x) - (a + b) * log1p(_ν1 * _x / _ν2) - logbeta(a, b)
    end
    fdistpdf_ref(ν1::Real, ν2::Real, x::Real) = exp(fdistlogpdf_ref(ν1, ν2, x))

    # the density has to stay accurate when either degree of freedom is huge: for large
    # `ν2` the term `ν1 * x / ν2` is tiny and must not be absorbed by an explicit `1 +`,
    # and for large `ν1` no two terms growing like `ν1 * log(ν1)` may cancel
    @testset "fdistlogpdf with ν$i = $ν" for (i, ν) in
        Iterators.product(1:2, (1.0e9, 1.0e12, 1.0e15, 1.0e18, 1.0e20))

        ν1, ν2 = i == 1 ? (ν, 5.0) : (1.0, ν)
        @testset "x = $x" for x in (0.5, 2.0, 6.6)
            @test fdistlogpdf(ν1, ν2, x) ≈ fdistlogpdf_ref(ν1, ν2, x) rtol = 1.0e-13
            @test fdistpdf(ν1, ν2, x) ≈ fdistpdf_ref(ν1, ν2, x) rtol = 1.0e-13
        end
    end
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

    # at zero the density behaves like `x^(ν1 / 2 - 1)`
    @testset "x = 0" begin
        @test fdistlogpdf(1, 1, 0.0) == Inf
        @test fdistlogpdf(2, 1, 0.0) == 0
        @test fdistlogpdf(10, 3, 0.0) == -Inf
    end

    # `NaN` propagates through the `xlogy` in the `x == 0` branch
    @testset "x = NaN" begin
        @test isnan(fdistlogpdf(3, 7, NaN))
        @test isnan(fdistpdf(3, 7, NaN))
    end
end
