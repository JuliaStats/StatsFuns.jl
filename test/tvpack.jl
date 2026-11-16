@testitem "TVPack" begin
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
        @test bvncdf(0.0, 100000.0, 0.0) ≈ 0.5
        @test bvncdf(100000.0, 100000.0, 0.0) ≈ 1.0
        @test bvncdf(-100000.0, 100000.0, 0.0) ≈ 0.0

        for x in -100.0:100.0, y in -100.0:100.0, r in -1.0:0.1:1.0
            @test !isnan(bvncdf(x, y, r))
        end
    end

    @testset "tvtcdf: trivariate t cdf" begin
        # r = [r12, r13, r23]; nu = 0 is the normal case
        h = [0.5, 0.2, -0.1]
        for nu in (0, 1, 3, 30), r in ([0.3, 0.2, 0.1], [-0.5, 0.4, 0.2], [0.9, 0.8, 0.7])
            # sum(abs, h) just above eps() skips the h ≈ 0 shortcut and exercises the numerical integration
            @test StatsFuns.tvtcdf(nu, fill(1.0e-15, 3), r) ≈ 1 / 8 + sum(asin, r) / 4π rtol = 1.0e-13

            # A huge h3 marginalizes out the third variable
            @test StatsFuns.tvtcdf(nu, [h[1], h[2], 1.0e15], r) ≈ StatsFuns.bvtcdf(nu, h[1], h[2], r[1]) rtol = 1.0e-13

            # Relabeling the variables does not change the probability
            p = StatsFuns.tvtcdf(nu, h, r)
            @test StatsFuns.tvtcdf(nu, h[[2, 1, 3]], r[[1, 3, 2]]) ≈ p rtol = 1.0e-14
            @test StatsFuns.tvtcdf(nu, h[[3, 2, 1]], r[[3, 2, 1]]) ≈ p rtol = 1.0e-14
        end
    end
end
