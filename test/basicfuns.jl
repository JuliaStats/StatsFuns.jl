using StatsFuns
using Base.Test
using Compat

# xlogx & xlogy

println("\ttesting xlogx & xlogy ...")

@test xlogx(0) === 0.0
@test xlogx(2) ≈ 2.0 * log(2.0)

@test xlogy(0, 1) === 0.0
@test xlogy(2, 3) ≈ 2.0 * log(3.0)

# logistic & logit

println("\ttesting logistic & logit ...")

@test logistic(2)        ≈ 1.0 / (1.0 + exp(-2.0))
@test logit(0.5)         ≈ 0.0
@test logit(logistic(2)) ≈ 2.0

# log1psq

println("\ttesting log1psq ...")

@test log1psq(0.0) ≈ 0.0
@test log1psq(1.0) ≈ log1p(1.0)
@test log1psq(2.0) ≈ log1p(4.0)

# log1pexp, log1mexp, log2mexp & logexpm1

println("\ttesting log1pexp ...")

@test log1pexp(2.0)    ≈ log(1.0 + exp(2.0))
@test log1pexp(-2.0)   ≈ log(1.0 + exp(-2.0))
@test log1pexp(10000)  ≈ 10000.0
@test log1pexp(-10000) ≈ 0.0

println("\ttesting log1mexp ...")

@test log1mexp(-1.0)  ≈ log1p(- exp(-1.0))
@test log1mexp(-10.0) ≈ log1p(- exp(-10.0))

println("\ttesting log2mexp ...")

@test log2mexp(0.0)  ≈ 0.0
@test log2mexp(-1.0) ≈ log(2.0 - exp(-1.0))

println("\ttesting logexpm1 ...")

@test logexpm1(2.0)            ≈  log(exp(2.0) - 1.0)
@test logexpm1(log1pexp(2.0))  ≈  2.0
@test logexpm1(log1pexp(-2.0)) ≈ -2.0

# log1pmx

println("\ttesting log1pmx ...")

@test log1pmx(0.0) ≈ 0.0
@test log1pmx(1.0) ≈ log(2.0) - 1.0
@test log1pmx(2.0) ≈ log(3.0) - 2.0

println("\ttesting logmxp1 ...")

@test logmxp1(1.0) ≈ 0.0
@test logmxp1(2.0) ≈ log(2.0) - 1.0
@test logmxp1(3.0) ≈ log(3.0) - 2.0

# logsumexp

println("\ttesting logsumexp ...")

@test logsumexp(2.0, 3.0)     ≈ log(exp(2.0) + exp(3.0))
@test logsumexp(10002, 10003) ≈ 10000 + logsumexp(2.0, 3.0)

@test logsumexp([1.0, 2.0, 3.0])          ≈ 3.40760596444438
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
        @test logsumexp(arguments) ≡ result
        @test logsumexp(arguments...) ≡ result
    end
end

# softmax

println("\ttesting softmax ...")

x = [1.0, 2.0, 3.0]
# Explicit conversion to Vector{Float64} can be romved when we drop 0.4 support
r = Vector{Float64}(@compat exp.(x) ./ sum(exp.(x)))
@test softmax([1.0, 2.0, 3.0]) ≈ r
softmax!(x)
@test x ≈ r
