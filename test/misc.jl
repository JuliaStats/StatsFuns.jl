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
