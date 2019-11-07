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
