using StatsFuns, Test

@testset "bvncdf: bivariate normal cdf" begin

    @test bvncdf(0.0, 0.0, 0.0) == 0.25
    @test bvncdf(0.0, Inf, 0.0) == 0.5
    @test bvncdf(Inf, 0.0, 0.0) == 0.5
    @test bvncdf(Inf, Inf, 0.0) == 1.0

    for r in -1.0:0.25:1.0
        for x in -10.0:0.5:10.0
            @test bvncdf(x, Inf, r) == normcdf(x)
            @test bvncdf(Inf, x, r) == normcdf(x)
            @test bvncdf(x, -Inf, r) == 0.0
            @test bvncdf(-Inf, x, r) == 0.0
        end
    end

    @test_throws DomainError bvncdf(0.0, 0.0, -2.0)
    @test_throws DomainError bvncdf(0.0, 0.0, 2.0)

    @test bvncdf(0.0, -100000.0, 0.0) ≈ 0.0
    @test bvncdf(0.0,  100000.0, 0.0) ≈ 0.5
    @test bvncdf(100000.0,  100000.0, 0.0000) ≈ 1.0
    @test bvncdf(-100000.0,  100000.0, 0.0000) ≈ 0.0

    for x in -100.0:100.0, y in -100.0:100.0, r in -1.0:0.1:1.0
        @test !isnan(bvncdf(x,y,r))
    end
end
