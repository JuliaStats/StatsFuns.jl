using SpecialFunctions, StatsFuns

@testset "logmvgamma" begin
    @testset "type behavior" for eltya in (Float32, Float64)
        p = rand(1:50)
        a = rand(eltya)
        @test typeof(logmvgamma(p, a)) == eltya
    end

    @testset "consistent with lgamma" for eltya in (Float32, Float64)
        #  Γ₁(a) = Γ(a), Γ₂(a) = √π Γ(a) Γ(a - 0.5), etc
        a = rand(eltya)
        @test logmvgamma(1, a) ≈ lgamma(a)
        @test logmvgamma(2, a) ≈    eltya(0.5logπ) + lgamma(a) + lgamma(a - eltya(0.5))
        @test logmvgamma(3, a) ≈ eltya((3/2)*logπ) + lgamma(a) + lgamma(a - eltya(0.5)) + lgamma(a - one(a))
    end

    @testset "consistent with itself" for eltya in (Float32, Float64)
        #  Γᵢ(a) = (π^{i-1/2}) Γ(a) Γᵢ₋₁(a - 0.5)
        for p in 1:50
            a = rand(eltya)
            @test logmvgamma(p, a) ≈ eltya((p/2-1/2)*logπ) + lgamma(a) + logmvgamma(p - 1, a - eltya(0.5))
        end
    end
end

@testset "logmvbeta" begin
    @testset "symmetry" for eltya in (Float32, Float64)
                            for eltyb in (Float32, Float64)
            #  Bᵢ(a, b) = Bᵢ(b, a)
            for p in 1:50
                a = rand(eltya)
                b = rand(eltyb)
                @test logmvbeta(p, a, b) ≈ logmvbeta(p, b, a)
            end
        end
    end

    @testset "consistent with logbeta" for eltya in (Float32, Float64)
                                           for eltyb in (Float32, Float64)
            #  B₁(a, b) = B(a, b)
            a = rand(eltya)
            b = rand(eltyb)
            @test logmvbeta(1, a, b) ≈ lgamma(a) + lgamma(b) - lgamma(a + b)
            #  Change to logbeta calls after next SpecialFunctions release
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
