@testitem "ForwardDiff gamma tails" tags = [:autodiff] begin
    using StatsFuns
    using ForwardDiff: ForwardDiff, Dual
    using Test

    # Closed forms of the Erlang tails, i.e., of the gamma tails with an integer shape `k`
    function erlang_logccdf(k::Integer, θ, x)
        u = x / θ
        return -u + log(sum(u^i / factorial(i) for i in 0:(k - 1)))
    end
    erlang_ccdf(k::Integer, θ, x) = exp(erlang_logccdf(k, θ, x))
    erlang_cdf(k::Integer, θ, x) = -expm1(erlang_logccdf(k, θ, x))
    erlang_logcdf(k::Integer, θ, x) = log1mexp(erlang_logccdf(k, θ, x))

    tails = (
        (gammacdf, erlang_cdf),
        (gammaccdf, erlang_ccdf),
        (gammalogcdf, erlang_logcdf),
        (gammalogccdf, erlang_logccdf),
    )

    erlang_cases = [
        (f, erlang, k, θ, x)
            for (f, erlang) in tails
            for k in (1, 2, 5)
            for θ in (1.0, 3.0)
            for x in (0.5, 1.0, 4.0, 12.0)
    ]
    @testset "Erlang closed form: $f($k, $θ, $x)" for (f, erlang, k, θ, x) in erlang_cases
        @test ForwardDiff.value(f(k, θ, Dual(x, 1.0))) ≈ f(k, θ, x)
        @test ForwardDiff.derivative(y -> f(k, θ, y), x) ≈
            ForwardDiff.derivative(y -> erlang(k, θ, y), x)
        @test ForwardDiff.derivative(t -> f(k, t, x), θ) ≈
            ForwardDiff.derivative(t -> erlang(k, t, x), θ)
        # `θ` and `x` with the same tag
        @test ForwardDiff.gradient(v -> f(k, v[1], v[2]), [θ, x]) ≈
            ForwardDiff.gradient(v -> erlang(k, v[1], v[2]), [θ, x])
        @test ForwardDiff.hessian(v -> f(k, v[1], v[2]), [θ, x]) ≈
            ForwardDiff.hessian(v -> erlang(k, v[1], v[2]), [θ, x])
        # `θ` and `x` with different tags
        @test ForwardDiff.derivative(t -> ForwardDiff.derivative(y -> f(k, t, y), x), θ) ≈
            ForwardDiff.derivative(t -> ForwardDiff.derivative(y -> erlang(k, t, y), x), θ)
        @test ForwardDiff.derivative(y -> ForwardDiff.derivative(t -> f(k, t, y), θ), x) ≈
            ForwardDiff.derivative(y -> ForwardDiff.derivative(t -> erlang(k, t, y), θ), x)
    end

    central_difference(g, t) = (g(t + 1.0e-6) - g(t - 1.0e-6)) / 2.0e-6

    difference_cases = [
        (f, k, θ, x)
            for (f, _) in tails
            for k in (0.5, 2.5)
            for θ in (1.0, 3.0)
            for x in (0.5, 1.0, 4.0)
    ]
    @testset "finite differences: $f($k, $θ, $x)" for (f, k, θ, x) in difference_cases
        @test ForwardDiff.derivative(y -> f(k, θ, y), x) ≈
            central_difference(y -> f(k, θ, y), x) rtol = 1.0e-6
        @test ForwardDiff.derivative(t -> f(k, t, x), θ) ≈
            central_difference(t -> f(k, t, x), θ) rtol = 1.0e-6
    end

    @testset "log tails where the tail probability underflows" begin
        # `gammaccdf(3, 1, 500)` underflows
        @test ForwardDiff.derivative(y -> gammalogccdf(3, 1.0, y), 500.0) ≈
            ForwardDiff.derivative(y -> erlang_logccdf(3, 1.0, y), 500.0)
        @test ForwardDiff.derivative(t -> gammalogccdf(3, t, 500.0), 1.0) ≈
            ForwardDiff.derivative(t -> erlang_logccdf(3, t, 500.0), 1.0)
        # `gammacdf(3, 1, x) ≈ x^3 / 6` underflows for `x = 1e-120`
        @test ForwardDiff.derivative(y -> gammalogcdf(3, 1.0, y), 1.0e-120) ≈ 3.0e120
        @test ForwardDiff.derivative(t -> gammalogcdf(3, t, 1.0e-120), 1.0) ≈ -3
    end

    boundary_cases = [(f, k) for (f, _) in tails for k in (0.5, 1.0, 3.0)]
    @testset "boundary of the support: $f($k, ...)" for (f, k) in boundary_cases
        @test iszero(ForwardDiff.gradient(v -> f(k, v[1], v[2]), [1.0, -1.0]))
        @test iszero(ForwardDiff.derivative(t -> f(k, t, 0.0), 1.0))
    end
    @testset "x-derivative at x = 0" begin
        @test iszero(ForwardDiff.derivative(y -> gammacdf(3.0, 1.0, y), 0.0))
        @test iszero(ForwardDiff.derivative(y -> gammaccdf(3.0, 1.0, y), 0.0))
        @test ForwardDiff.derivative(y -> gammacdf(1.0, 2.0, y), 0.0) ≈ 0.5
        @test ForwardDiff.derivative(y -> gammaccdf(1.0, 2.0, y), 0.0) ≈ -0.5
    end

    # The Poisson tails promote `λ` and the shape to `Dual`s, which `_gammalogcdf` and
    # `_gammalogccdf` do not accept (JuliaStats/StatsFuns.jl#161).
    poisson_cases = [(λ, n) for λ in (0.5, 1.3, 7.0) for n in (0, 3, 10)]
    @testset "Poisson tails differentiate (λ=$λ, n=$n)" for (λ, n) in poisson_cases
        @test ForwardDiff.derivative(l -> poiscdf(l, n), λ) ≈ -poispdf(λ, n)
        @test ForwardDiff.derivative(l -> poisccdf(l, n), λ) ≈ poispdf(λ, n)
        @test ForwardDiff.derivative(l -> poislogcdf(l, n), λ) ≈
            -poispdf(λ, n) / poiscdf(λ, n)
        @test ForwardDiff.derivative(l -> poislogccdf(l, n), λ) ≈
            poispdf(λ, n) / poisccdf(λ, n)
    end

    @testset "$f with a shape tangent throws" for (f, _) in tails
        @test_throws ArgumentError ForwardDiff.derivative(k -> f(k, 1.0, 2.0), 3.0)
        @test_throws ArgumentError ForwardDiff.gradient(v -> f(v[1], v[2], 2.0), [3.0, 1.0])
        @test_throws ArgumentError ForwardDiff.gradient(v -> f(v[1], 1.0, v[2]), [3.0, 2.0])
        @test_throws ArgumentError ForwardDiff.gradient(
            v -> f(v[1], v[2], v[3]),
            [3.0, 1.0, 2.0],
        )
    end

    @testset "$f with a promoted shape (zero partials)" for (f, erlang) in tails
        # e.g. `Distributions.Gamma(2.0, θ)` promotes the shape to the type of `θ`
        dθ = ForwardDiff.derivative(t -> erlang(2, t, 1.5), 1.0)
        dx = ForwardDiff.derivative(y -> erlang(2, 1.0, y), 1.5)
        @test ForwardDiff.derivative(t -> f(oftype(t, 2), t, 1.5), 1.0) ≈ dθ
        @test ForwardDiff.derivative(y -> f(oftype(y, 2), 1.0, y), 1.5) ≈ dx
        @test ForwardDiff.gradient(v -> f(2 + 0 * v[1], v[1], v[2]), [1.0, 1.5]) ≈ [dθ, dx]
        @test ForwardDiff.partials(f(Dual(2.0, 0.0), 1.0, Dual(1.5, 1.0)))[1] ≈ dx
        # a zero-partial shape with a different tag than `θ`
        @test ForwardDiff.derivative(
            t -> ForwardDiff.derivative(s -> f(2 + 0 * s, t, 1.5), 1.0),
            1.0,
        ) ==
            0
        @test ForwardDiff.derivative(
            s -> ForwardDiff.derivative(t -> f(2 + 0 * s, t, 1.5), 1.0),
            1.0,
        ) ==
            0
    end
end
