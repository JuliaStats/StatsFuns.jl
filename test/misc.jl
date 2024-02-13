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
        @test logmvgamma(2, a) ≈    eltya(0.5logπ) + loggamma(a) + loggamma(a - eltya(0.5))
        @test logmvgamma(3, a) ≈ eltya((3/2)*logπ) + loggamma(a) + loggamma(a - eltya(0.5)) + loggamma(a - one(a))
    end

    @testset "consistent with itself" for eltya in (Float32, Float64)
        #  Γᵢ(a) = (π^{i-1/2}) Γ(a) Γᵢ₋₁(a - 0.5)
        for p in 1:50
            a = rand(eltya) + p # add p since loggamma is only define for positive arguments
            @test logmvgamma(p, a) ≈ eltya((p/2-1/2)*logπ) + loggamma(a) + logmvgamma(p - 1, a - eltya(0.5))
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
            # Older SpecialFunctions + Julia versions require slightly larger tolerance
            if VERSION < v"1.3" && promote_type(eltya, eltyb) === Float64
                @test logmvbeta(1, a, b) ≈ logbeta(a, b) rtol = 10 * sqrt(eps(Float64))
            else
                @test logmvbeta(1, a, b) ≈ logbeta(a, b)
            end
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
    @test gammalogcdf(42648.50647826826, 2.2498007956420723e-5, 0.6991377135675367) ≈ -1933.2698959040617410
end

# https://github.com/JuliaStats/StatsFuns.jl/issues/154
@testset "tvdistinvcdf: numerical issue" begin
    @test isnan(@inferred(tdistinvcdf(0, 0.975)))
end
