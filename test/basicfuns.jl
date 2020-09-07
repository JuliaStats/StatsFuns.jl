using StatsFuns, Test

@testset "xlogx & xlogy" begin
    @test iszero(xlogx(0))
    @test xlogx(2) ≈ 2.0 * log(2.0)
    @test_throws DomainError xlogx(-1)
    @test isnan(xlogx(NaN))

    @test iszero(xlogy(0, 1))
    @test xlogy(2, 3) ≈ 2.0 * log(3.0)
    @test_throws DomainError xlogy(1, -1)
    @test isnan(xlogy(NaN, 2))
    @test isnan(xlogy(2, NaN))
    @test isnan(xlogy(0, NaN))

    # Since we allow complex/negative values, test for them. See comments in:
    # https://github.com/JuliaStats/StatsFuns.jl/pull/95

    @test xlogx(1 + im) == (1 + im) * log(1 + im)
    @test isnan(xlogx(NaN + im))
    @test isnan(xlogx(1 + NaN * im))

    @test xlogy(-2, 3) == -xlogy(2, 3)
    @test xlogy(1 + im, 3) == (1 + im) * log(3)
    @test xlogy(1 + im, 2 + im) == (1 + im) * log(2 + im)
    @test isnan(xlogy(1 + NaN * im, -1 + im))
    @test isnan(xlogy(0, -1 + NaN * im))
    @test isnan(xlogy(Inf + im * NaN, 1))
    @test isnan(xlogy(0 + im * 0, NaN))
    @test iszero(xlogy(0 + im * 0, 0 + im * Inf))
end

@testset "logistic & logit" begin
    @test logistic(2) ≈ 1.0 / (1.0 + exp(-2.0))
    @test logistic(-750.0) === 0.0
    @test logistic(-740.0) > 0.0
    @test logistic(+36.0) < 1.0
    @test logistic(+750.0) === 1.0
    @test iszero(logit(0.5))
    @test logit(logistic(2)) ≈ 2.0
end

@testset "log1psq" begin
    @test iszero(log1psq(0.0))
    @test log1psq(1.0) ≈ log1p(1.0)
    @test log1psq(2.0) ≈ log1p(4.0)
end

# log1pexp, log1mexp, log2mexp & logexpm1

@testset "log1pexp" begin
    @test log1pexp(2.0)    ≈ log(1.0 + exp(2.0))
    @test log1pexp(-2.0)   ≈ log(1.0 + exp(-2.0))
    @test log1pexp(10000)  ≈ 10000.0
    @test log1pexp(-10000) ≈ 0.0

    @test log1pexp(2f0)      ≈ log(1f0 + exp(2f0))
    @test log1pexp(-2f0)     ≈ log(1f0 + exp(-2f0))
    @test log1pexp(10000f0)  ≈ 10000f0
    @test log1pexp(-10000f0) ≈ 0f0
end

@testset "log1mexp" begin
    @test log1mexp(-1.0)  ≈ log1p(- exp(-1.0))
    @test log1mexp(-10.0) ≈ log1p(- exp(-10.0))
end

@testset "log2mexp" begin
    @test log2mexp(0.0)  ≈ 0.0
    @test log2mexp(-1.0) ≈ log(2.0 - exp(-1.0))
end

@testset "logexpm1" begin
    @test logexpm1(2.0)            ≈  log(exp(2.0) - 1.0)
    @test logexpm1(log1pexp(2.0))  ≈  2.0
    @test logexpm1(log1pexp(-2.0)) ≈ -2.0

    @test logexpm1(2f0)            ≈  log(exp(2f0) - 1f0)
    @test logexpm1(log1pexp(2f0))  ≈  2f0
    @test logexpm1(log1pexp(-2f0)) ≈ -2f0
end

@testset "log1pmx" begin
    @test iszero(log1pmx(0.0))
    @test log1pmx(1.0) ≈ log(2.0) - 1.0
    @test log1pmx(2.0) ≈ log(3.0) - 2.0
end

@testset "logmxp1" begin
    @test iszero(logmxp1(1.0))
    @test logmxp1(2.0) ≈ log(2.0) - 1.0
    @test logmxp1(3.0) ≈ log(3.0) - 2.0
end

@testset "logsumexp" begin
    @test logaddexp(2.0, 3.0)     ≈ log(exp(2.0) + exp(3.0))
    @test logaddexp(10002, 10003) ≈ 10000 + logaddexp(2.0, 3.0)

    @test @inferred(logsumexp([1.0])) == 1.0
    @test @inferred(logsumexp((x for x in [1.0]))) == 1.0
    @test @inferred(logsumexp([1.0, 2.0, 3.0])) ≈ 3.40760596444438
    @test @inferred(logsumexp((1.0, 2.0, 3.0))) ≈ 3.40760596444438
    @test logsumexp([1.0, 2.0, 3.0] .+ 1000.) ≈ 1003.40760596444438

    @test @inferred(logsumexp([[1.0, 2.0, 3.0] [1.0, 2.0, 3.0] .+ 1000.]; dims=1)) ≈ [3.40760596444438 1003.40760596444438]
    @test @inferred(logsumexp([[1.0 2.0 3.0]; [1.0 2.0 3.0] .+ 1000.]; dims=2)) ≈ [3.40760596444438, 1003.40760596444438]
    @test @inferred(logsumexp([[1.0, 2.0, 3.0] [1.0, 2.0, 3.0] .+ 1000.]; dims=[1,2])) ≈ [1003.4076059644444]

    # check underflow
    @test logsumexp([1e-20, log(1e-20)]) ≈ 2e-20

    let cases = [([-Inf, -Inf], -Inf),   # correct handling of all -Inf
                 ([-Inf, -Inf32], -Inf), # promotion
                 ([-Inf32, -Inf32], -Inf32), # Float32
                 ([-Inf, Inf], Inf),
                 ([-Inf, 9.0], 9.0),
                 ([Inf, 9.0], Inf),
                 ([0, 0], log(2.0))] # non-float arguments
        for (arguments, result) in cases
            @test logaddexp(arguments...) ≡ result
            @test logsumexp(arguments) ≡ result
        end
    end

    @test isnan(logsubexp(Inf, Inf))
    @test isnan(logsubexp(-Inf, -Inf))
    @test logsubexp(Inf, 9.0) ≡ Inf
    @test logsubexp(-Inf, 9.0) ≡ 9.0
    @test logsubexp(1f2, 1f2) ≡ -Inf32
    @test logsubexp(0, 0) ≡ -Inf
    @test logsubexp(3, 2) ≈ 2.541324854612918108978

    # NaN propagation
    @test isnan(logaddexp(NaN, 9.0))
    @test isnan(logaddexp(NaN, Inf))
    @test isnan(logaddexp(NaN, -Inf))

    @test isnan(logsubexp(NaN, 9.0))
    @test isnan(logsubexp(NaN, Inf))
    @test isnan(logsubexp(NaN, -Inf))

    @test isnan(logsumexp([NaN, 9.0]))
    @test isnan(logsumexp([NaN, Inf]))
    @test isnan(logsumexp([NaN, -Inf]))

    # logsumexp with general iterables (issue #63)
    xs = range(-500, stop = 10, length = 1000)
    @test @inferred(logsumexp(x for x in xs)) == logsumexp(xs)
end

@testset "softmax" begin
    x = [1.0, 2.0, 3.0]
    r = exp.(x) ./ sum(exp.(x))
    @test softmax([1.0, 2.0, 3.0]) ≈ r
    softmax!(x)
    @test x ≈ r
end
