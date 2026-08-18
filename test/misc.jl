@testitem "Misc" begin
    using SpecialFunctions, StatsFuns
    using Test

    @testset "logmvgamma" begin
        @testset "type behavior" for eltya in (Float32, Float64)
            p = rand(1:50)
            a = rand(eltya)
            # add p since loggamma is only define for positive arguments
            @test typeof(logmvgamma(p, a + p)) == eltya
        end

        @testset "consistent with loggamma" for eltya in (Float32, Float64)
            #  Γ₁(a) = Γ(a), Γ₂(a) = √π Γ(a) Γ(a - 0.5), etc
            a = rand(eltya) + 1 # add one since loggamma is only define for positive arguments
            @test logmvgamma(1, a) ≈ loggamma(a)
            @test logmvgamma(2, a) ≈ eltya(0.5logπ) + loggamma(a) + loggamma(a - eltya(0.5))
            @test logmvgamma(3, a) ≈ eltya((3 / 2) * logπ) + loggamma(a) + loggamma(a - eltya(0.5)) + loggamma(a - one(a))
        end

        @testset "consistent with itself" for eltya in (Float32, Float64)
            #  Γᵢ(a) = (π^{i-1/2}) Γ(a) Γᵢ₋₁(a - 0.5)
            for p in 1:50
                a = rand(eltya) + p # add p since loggamma is only define for positive arguments
                @test logmvgamma(p, a) ≈ eltya((p / 2 - 1 / 2) * logπ) + loggamma(a) + logmvgamma(p - 1, a - eltya(0.5))
            end
        end
    end

    @testset "logmvbeta" begin
        @testset "symmetry" for eltya in (Float32, Float64)
            for eltyb in (Float32, Float64)
                #  Bᵢ(a, b) = Bᵢ(b, a)
                for p in 1:50
                    a = rand(eltya) + p
                    b = rand(eltyb) + p
                    @test logmvbeta(p, a, b) ≈ logmvbeta(p, b, a)
                end
            end
        end

        @testset "consistent with logbeta" for eltya in (Float32, Float64)
            for eltyb in (Float32, Float64)
                #  B₁(a, b) = B(a, b)
                a = rand(eltya)
                b = rand(eltyb)
                @test logmvbeta(1, a, b) ≈ logbeta(a, b)
            end
        end

        @testset "type promotion behaves" for eltya in (Float32, Float64)
            for eltyb in (Float32, Float64)
                a = rand(eltya)
                b = rand(eltyb)
                T = Base.promote_type(eltya, eltyb)
                @test typeof(logmvbeta(1, a, b)) == T
            end
        end
    end

    # https://github.com/JuliaStats/StatsFuns.jl/issues/115
    @testset "support of binomial distribution" begin
        @test iszero(binompdf(1, 0.5, prevfloat(1.0)))
        @test iszero(binompdf(1, 0.5, nextfloat(1.0)))
        @test binomlogpdf(1, 0.5, prevfloat(1.0)) == -Inf
        @test binomlogpdf(1, 0.5, nextfloat(1.0)) == -Inf
    end

    @testset "binom special cases" begin
        for (n, p, k) in ((5, 0.0, 0), (5, 1.0, 5))
            @test iszero(binomlogpdf(n, p, k))
            @test isone(binompdf(n, p, k))
        end
    end

    @testset "lstirling_asym" begin
        # can test for equality here because the lhs is the way the value is created
        @test Float32(lstirling_asym(1.0)) == @inferred lstirling_asym(1.0f0)
        # for 64.0f0 the expansion is used but for 64.0 the BigFloat value is rounded
        @test Float32(lstirling_asym(64.0)) ≈ @inferred lstirling_asym(64.0f0)
    end

    # https://github.com/JuliaStats/StatsFuns.jl/issues/143
    # https://github.com/JuliaMath/HypergeometricFunctions.jl/issues/47
    @testset "betalogcdf: numerical issue" begin
        # Mathematica: N[Log[CDF[BetaDistribution[6041, 2496], 1/10]], 10]
        @test betalogcdf(6041, 2496, 0.1) ≈ -9020.029401
        @test betainvlogcdf(6041, 2496, betalogcdf(6041, 2496, 0.1)) ≈ 0.1
    end

    # https://github.com/JuliaStats/StatsFuns.jl/issues/150
    @testset "gammalogcdf: numerical issue" begin
        @test gammalogcdf(42648.50647826826, 2.2498007956420723e-5, 0.6991377135675367) ≈ -1933.269895904061741
    end

    # https://github.com/JuliaStats/StatsFuns.jl/issues/154
    @testset "tvdistinvcdf: numerical issue" begin
        @test isnan(@inferred(tdistinvcdf(0, 0.975)))
    end
end

@testitem "tdistinvcdf" begin
    using StatsFuns
    using Test

    @testset "reference values" begin
        # reference values computed with a 512-bit MPFR Newton iteration on the
        # regularized incomplete beta representation of the cdf
        for (ν, p, t) in (
                (0.35, 0.499, -0.0041388393417127554),
                (0.5, 0.3, -1.0095258786071661),
                (0.5, 1.0e-8, -1.02849115631634e15),
                (1.0, 1.0e-100, -3.1830988618379064e99),
                (2.5, 1.0e-12, -55306.174076515817),
                (5.0, 1.0e-8, -62.40450611096729),
                (5.0, 0.95, 2.0150483733330233),
                (12.0, 0.3, -0.53861766820191637),
                (50.0, 0.999, 3.261409055798318),
                (300.0, 0.025, -1.9679030112610869),
                (1000.0, 0.975, 1.9623390808264081),
                (1.0e4, 1.0e-20, -9.282474153254304),
                (1.0e6, 1 - 1.0e-12, 7.0345756932732169),
            )
            @test tdistinvcdf(ν, p) ≈ t rtol = 1.0e-13
            @test tdistinvccdf(ν, p) == -tdistinvcdf(ν, p)
        end
    end

    @testset "round trips" begin
        # t -> cdf -> invcdf on a grid where the cdf itself is accurate; deep
        # tails are covered by the reference values above since tdistcdf
        # currently loses precision there
        # the quantile is always represented through its own tail: going through
        # the complementary probability is ill-conditioned for any implementation
        for ν in (0.5, 1.0, 2.5, 5.0, 20.0, 100.0, 1.0e3, 1.0e6),
                t in (-8.0, -3.0, -0.5, 0.0, 1.0, 6.0)

            if t <= 0
                p = tdistcdf(ν, t)
                0 < p < 1 || continue
                @test tdistinvcdf(ν, p) ≈ t atol = 1.0e-14 rtol = 1.0e-11
            else
                q = tdistccdf(ν, t)
                0 < q < 1 || continue
                @test tdistinvccdf(ν, q) ≈ t atol = 1.0e-14 rtol = 1.0e-11
            end
        end
    end

    @testset "edge cases and types" begin
        @test tdistinvcdf(5, 0.0) == -Inf
        @test tdistinvcdf(5, 1.0) == Inf
        @test tdistinvcdf(5, 0.5) === 0.0
        @test isnan(tdistinvcdf(5.0, NaN))
        @test isnan(tdistinvcdf(NaN, 0.5))
        @test isnan(tdistinvcdf(5.0, -0.1))
        @test isnan(tdistinvcdf(5.0, 1.1))
        @test isnan(tdistinvcdf(-1.0, 0.5))
        @test tdistinvcdf(Inf, 0.975) == norminvcdf(0.975)
        @test @inferred(tdistinvcdf(1, 0.75f0)) isa Float32
        @test @inferred(tdistinvcdf(Float16(1), Float16(0.75))) isa Float16
        @test @inferred(tdistinvcdf(1, 0.75)) isa Float64
        # same-type non-IEEE arguments use the same kernel
        @test tdistinvcdf(1 // 2, 1 // 100) == tdistinvcdf(0.5, 0.01)
        @test @inferred(tdistinvcdf(5, 1)) == Inf
        # StatsFuns#228: deep tails for small ν used to underflow to -Inf
        @test isfinite(tdistinvcdf(0.5, 1.0e-8))
        # symmetry (1 - 0.75 is exact in binary)
        @test tdistinvcdf(3.5, 0.25) == -tdistinvcdf(3.5, 0.75)
    end
end
