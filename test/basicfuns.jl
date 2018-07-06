using StatsFuns, Test

@testset "xlogx & xlogy" begin
    @test iszero(xlogx(0))
    @test xlogx(2) ≈ 2.0 * log(2.0)

    @test iszero(xlogy(0, 1))
    @test xlogy(2, 3) ≈ 2.0 * log(3.0)
end

@testset "logistic & logit" begin
    @test logistic(2)        ≈ 1.0 / (1.0 + exp(-2.0))
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

    @test logsumexp([1.0, 2.0, 3.0])          ≈ 3.40760596444438
    @test logsumexp((1.0, 2.0, 3.0))          ≈ 3.40760596444438
    @test logsumexp([1.0, 2.0, 3.0] .+ 1000.) ≈ 1003.40760596444438

    let cases = [([-Inf, -Inf], -Inf),   # correct handling of all -Inf
                 ([-Inf, -Inf32], -Inf), # promotion
                 ([-Inf32, -Inf32], -Inf32), # Float32
                 ([-Inf, Inf], Inf),
                 ([-Inf, 9.0], 9.0),
                 ([Inf, 9.0], Inf),
                 ([NaN, 9.0], NaN),  # NaN propagation
                 ([NaN, Inf], NaN),  # NaN propagation
                 ([NaN, -Inf], NaN), # NaN propagation
                 ([0, 0], log(2.0))] # non-float arguments
        for (arguments, result) in cases
            @test logaddexp(arguments...) ≡ result
            @test logsumexp(arguments) ≡ result
        end
    end
end

@testset "softmax" begin
    x = [1.0, 2.0, 3.0]
    r = exp.(x) ./ sum(exp.(x))
    @test softmax([1.0, 2.0, 3.0]) ≈ r
    softmax!(x)
    @test x ≈ r
end
