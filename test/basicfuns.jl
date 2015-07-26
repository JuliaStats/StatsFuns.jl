using StatsFuns
using Base.Test

# xlogx & xlogy

@test xlogx(0) === 0.0
@test_approx_eq xlogx(2) 2.0 * log(2.0)

@test xlogy(0, 1) === 0.0
@test_approx_eq xlogy(2, 3) 2.0 * log(3.0)

# logistic & logit

@test_approx_eq logistic(2) 1.0 / (1.0 + exp(-2.0))
@test_approx_eq logit(0.5) 0.0
@test_approx_eq logit(logistic(2)) 2.0

# log1psq

@test_approx_eq log1psq(0.0) 0.0
@test_approx_eq log1psq(1.0) log1p(1.0)
@test_approx_eq log1psq(2.0) log1p(4.0)

# log1pexp, log1mexp, log2mexp & logexpm1

@test_approx_eq log1pexp(2.0) log(1.0 + exp(2.0))
@test_approx_eq log1pexp(-2.0) log(1.0 + exp(-2.0))
@test_approx_eq log1pexp(10000) 10000.0
@test_approx_eq log1pexp(-10000) 0.0

@test_approx_eq log1mexp(-1.0) log1p(- exp(-1.0))
@test_approx_eq log1mexp(-10.0) log1p(- exp(-10.0))

@test_approx_eq log2mexp(0.0) 0.0
@test_approx_eq log2mexp(-1.0) log(2.0 - exp(-1.0))

@test_approx_eq logexpm1(2.0) log(exp(2.0) - 1.0)
@test_approx_eq logexpm1(log1pexp(2.0)) 2.0
@test_approx_eq logexpm1(log1pexp(-2.0)) -2.0

# log1pmx

@test_approx_eq log1pmx(0.0) 0.0
@test_approx_eq log1pmx(1.0) log(2.0) - 1.0
@test_approx_eq log1pmx(2.0) log(3.0) - 2.0

@test_approx_eq logmxp1(1.0) 0.0
@test_approx_eq logmxp1(2.0) log(2.0) - 1.0
@test_approx_eq logmxp1(3.0) log(3.0) - 2.0

# logsumexp

@test_approx_eq logsumexp(2.0, 3.0) log(exp(2.0) + exp(3.0))
@test_approx_eq logsumexp(10002, 10003) 10000 + logsumexp(2.0, 3.0)

@test_approx_eq logsumexp([1.0, 2.0, 3.0]) 3.40760596444438
@test_approx_eq logsumexp([1.0, 2.0, 3.0] .+ 1000.) 1003.40760596444438

# softmax

x = [1.0, 2.0, 3.0]
r = exp(x) ./ sum(exp(x))
@test_approx_eq softmax([1.0, 2.0, 3.0]) r
softmax!(x)
@test_approx_eq x r
