@testitem "ForwardDiff gamma tails" tags = [:autodiff] begin
    using StatsFuns
    using ForwardDiff: ForwardDiff, Dual
    using Test

    configurations = [
        (k, θ, x)
        for k in (0.5, 1.0, 3.0, 9.0)
        for θ in (1.0, 3.0)
        for x in (0.05, 1.0, 4.0, 12.0)
    ]

    @testset "value survives the tangent (k=$k, θ=$θ, x=$x)" for (k, θ, x) in configurations
        @test gammacdf(k, θ, Dual(x, 1.0)).value ≈ gammacdf(k, θ, x)
        @test gammaccdf(k, θ, Dual(x, 1.0)).value ≈ gammaccdf(k, θ, x)
        @test gammalogcdf(k, θ, Dual(x, 1.0)).value ≈ gammalogcdf(k, θ, x)
        @test gammalogccdf(k, θ, Dual(x, 1.0)).value ≈ gammalogccdf(k, θ, x)
    end

    @testset "x-derivative is the density (k=$k, θ=$θ, x=$x)" for (
            k,
            θ,
            x,
        ) in configurations
        density = gammapdf(k, θ, x)
        @test ForwardDiff.derivative(y -> gammacdf(k, θ, y), x) ≈ density
        @test ForwardDiff.derivative(y -> gammaccdf(k, θ, y), x) ≈ -density
        @test ForwardDiff.derivative(y -> gammalogcdf(k, θ, y), x) ≈
            density / gammacdf(k, θ, x)
        @test ForwardDiff.derivative(y -> gammalogccdf(k, θ, y), x) ≈
            -density / gammaccdf(k, θ, x)
    end

    @testset "θ is a scale (k=$k, θ=$θ, x=$x)" for (k, θ, x) in configurations
        density = gammapdf(k, θ, x)
        @test ForwardDiff.derivative(t -> gammacdf(k, t, x), θ) ≈ -(x / θ) * density
        @test ForwardDiff.derivative(t -> gammalogccdf(k, t, x), θ) ≈
            (x / θ) * density / gammaccdf(k, θ, x)
    end

    # The Poisson tails reach `gamma_inc` through an integer shape, which is the
    # dispatch failure reported in JuliaStats/StatsFuns.jl#161.
    @testset "Poisson tails differentiate (λ=$λ, n=$n)" for λ in (0.5, 1.3, 7.0), n in (
            0,
            3,
            10,
        )

        @test ForwardDiff.derivative(l -> poiscdf(l, n), λ) ≈ -poispdf(λ, n)
        @test ForwardDiff.derivative(l -> poisccdf(l, n), λ) ≈ poispdf(λ, n)
        @test ForwardDiff.derivative(l -> poislogcdf(l, n), λ) ≈
            -poispdf(λ, n) / poiscdf(λ, n)
        @test ForwardDiff.derivative(l -> poislogccdf(l, n), λ) ≈
            poispdf(λ, n) / poisccdf(λ, n)
    end

    @testset "a shape tangent is rejected" begin
        @test_throws ArgumentError ForwardDiff.derivative(k -> gammacdf(k, 1.0, 2.0), 3.0)
        @test_throws ArgumentError ForwardDiff.derivative(
            k -> gammalogccdf(k, 1.0, 2.0),
            3.0,
        )
    end

    @testset "second derivatives nest" begin
        curvature = ForwardDiff.derivative(
            y -> ForwardDiff.derivative(z -> gammacdf(3.0, 1.0, z), y),
            4.0,
        )
        @test curvature ≈ ForwardDiff.derivative(y -> gammapdf(3.0, 1.0, y), 4.0)
    end
end
