using SpecialFunctions, StatsFuns

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
            @test logmvbeta(1, a, b) ≈ logbeta(a, b)
        end
    end

    @testset "type promotion behaves" for eltya in (Float32, Float64)
                                      for eltyb in (Float32, Float64)
            a = rand(eltya)
            b = rand(eltyb)
            T = Base.promote_eltype(eltya, eltyb)
            @test typeof(logmvbeta(1, a, b)) == T
        end
    end
end

@testset "binomlogpdf" begin
    @testset "identity" begin

        for p in [rand(50); 0; 1]
            for n in 0:50
                for k in Any[-50:50; randn(10)]
                    @test binomlogpdf(n, Float32(p), k) ≈ log(binomlogpdf(n, p, k))
                end
            end
        end
    end
end
