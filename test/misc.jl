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

    # https://github.com/JuliaStats/StatsFuns.jl/issues/224
    @testset "fdist: boundary values" begin
        # outside the support, where the beta variate is not defined
        @test @inferred(fdistlogpdf(3, 7, -1.0)) == -Inf
        @test @inferred(fdistlogpdf(3, 7, -Inf)) == -Inf
        @test @inferred(fdistpdf(3, 7, -1.0)) == 0

        # the previous formulation of the density subtracted two infinities here
        @test @inferred(fdistlogpdf(3, 7, Inf)) == -Inf
        @test @inferred(fdistpdf(3, 7, Inf)) == 0
        @test @inferred(fdistcdf(3, 7, Inf)) == 1
        @test @inferred(fdistccdf(3, 7, Inf)) == 0

        # at zero the density behaves like `x^(ν1 / 2 - 1)`; `-0.0` is not `< 0` and takes the
        # same branch as `0.0`
        @test @inferred(fdistlogpdf(1, 1, 0.0)) == Inf
        @test @inferred(fdistlogpdf(2, 1, 0.0)) == 0
        @test @inferred(fdistlogpdf(10, 3, 0.0)) == -Inf
        @test @inferred(fdistlogpdf(10, 3, -0.0)) == -Inf
        @test @inferred(fdistcdf(3, 7, 0.0)) == 0

        @test isnan(@inferred(fdistlogpdf(3, 7, NaN)))
        @test isnan(@inferred(fdistpdf(3, 7, NaN)))
        @test isnan(@inferred(fdistcdf(3, 7, NaN)))
    end

    # `x < 0`, `x == 0` and `x > 0` take three different branches of `fdistlogpdf`, which all
    # have to return `float(promote_type(...))` for every mix of argument types
    @testset "fdist: type stability" begin
        @testset "(ν1, ν2, x) = ($ν1, $ν2, $x)" for (ν1, ν2) in
                ((3, 7), (3, 7.0f0), (3 // 1, 7), (3.0, Float16(7)), (big(3), 7)),
                x in (-1, 0, 2)

            T = float(Base.promote_typeof(ν1, ν2, x))
            @test @inferred(fdistlogpdf(ν1, ν2, x)) isa T
            @test @inferred(fdistpdf(ν1, ν2, x)) isa T
        end

        # `betacdf` and friends do not support `BigFloat`
        @testset "(ν1, ν2, x) = ($ν1, $ν2, $x)" for (ν1, ν2) in
                ((3, 7), (3, 7.0f0), (3 // 1, 7), (3.0, Float16(7))),
                x in (-1, 0, 2)

            T = float(Base.promote_typeof(ν1, ν2, x))
            @test @inferred(fdistcdf(ν1, ν2, x)) isa T
            @test @inferred(fdistccdf(ν1, ν2, x)) isa T
            @test @inferred(fdistlogcdf(ν1, ν2, x)) isa T
            @test @inferred(fdistlogccdf(ν1, ν2, x)) isa T
        end
    end
end

# https://github.com/JuliaStats/StatsFuns.jl/issues/224
@testitem "fdist: numerical issue" setup = [FDistRef] begin
    using StatsFuns
    using Test

    # the density has to stay accurate when either degree of freedom is huge: for large `ν2`
    # the term `ν1 * x / ν2` is tiny and must not be absorbed by an explicit `1 +`, and for
    # large `ν1` no two terms growing like `ν1 * log(ν1)` may cancel
    νs = (1.0e9, 1.0e12, 1.0e15, 1.0e18, 1.0e20)
    @testset "fdistlogpdf with (ν1, ν2) = ($ν1, $ν2)" for (ν1, ν2) in
        Iterators.flatten(((ν, 5.0), (1.0, ν)) for ν in νs)

        # `x = ν2 / ν1` is the ridge where the textbook form cancels, the two extremes are
        # where the ratio `ν1 * x / ν2` under- and overflows
        @testset "x = $x" for x in
            (1.0e-300, 1.0e-8, 0.5, 2.0, 6.6, 1.0e8, 1.0e300, 0.99 * ν2 / ν1, ν2 / ν1)

            @test @inferred(fdistlogpdf(ν1, ν2, x))::Float64 ≈ FDistRef.logpdf(ν1, ν2, x) rtol = 1.0e-12
        end
    end

    # with both degrees of freedom of `floatmax` scale the product `ν1 * x` overflows where the
    # ratio `(ν1 / ν2) * x` and the density itself are ordinary
    @testset "huge (ν1, ν2) = ($ν1, $ν2)" for (ν1, ν2) in ((1.0e200, 1.0e300), (1.0e300, 1.0e200))
        @testset "x = $x" for x in
            (1.0e-300, 1.0e-8, 0.5, 2.0, 1.0e8, 1.0e200, floatmax(Float64), ν2 / ν1)

            @test @inferred(fdistlogpdf(ν1, ν2, x))::Float64 ≈ FDistRef.logpdf(ν1, ν2, x) rtol = 1.0e-12
        end
    end

    # `ν2 / (ν1 * x)` overflows for `x < ν2 / (ν1 * floatmax(T))` - an ordinary argument in
    # `Float16` - and used to take the density and the cdf with it
    @testset "small x = $x" for x in (1.0e-300, 1.0f-38, Float16(1.0e-4))
        T = typeof(x)
        ν1, ν2 = one(T), T(7)
        @test @inferred(fdistlogpdf(ν1, ν2, x))::T ≈ FDistRef.logpdf(ν1, ν2, x) rtol = eps(T)^(3 // 4)

        ref = FDistRef.logcdf_smallx(ν1, ν2, x)
        @test @inferred(fdistlogcdf(ν1, ν2, x))::T ≈ ref rtol = eps(T)^(2 // 3)
        # `exp` turns the absolute error of the log-cdf into a relative one
        @test @inferred(fdistcdf(ν1, ν2, x))::T ≈ T(exp(ref)) rtol = eps(T)^(2 // 3) * max(1, abs(Float64(ref)))
    end

    # symmetrically, `ν1 * x` overflows for the largest `x`
    @testset "large x = $x" for x in (floatmax(Float64), floatmax(Float32), floatmax(Float16))
        T = typeof(x)
        ν1, ν2 = T(4), T(100)
        @test @inferred(fdistlogpdf(ν1, ν2, x))::T ≈ T(FDistRef.logpdf(ν1, ν2, x)) rtol = eps(T)^(3 // 4)
    end

    # the beta variate of the cdf overflows for a huge `ν2` if formed as `ν2 / (ν1 * x)`, and for a
    # huge `x` if formed as `(ν1 * x) / (ν1 * x + ν2)` - the latter takes `u` to `0`, not `1`
    @test fdistcdf(1.0, floatmax(Float64), floatmax(Float64)) ==
        betacdf(0.5, floatmax(Float64) / 2, 0.5)
    @test fdistcdf(1.0, 1.0e308, 1.0e308) == 1
    @test fdistcdf(1.0, 1.0e308, 0.9e308) == 1
end
