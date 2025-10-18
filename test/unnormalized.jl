using StatsFuns
using SpecialFunctions
using Test

@testset "logupdf + logulikelihood" begin
    @testset "optimized" begin
        # Beta distribution
        α = 0.1
        β = 2.5
        x = 0.7
        @test betalogupdf(α, β, x) ≈ betalogpdf(α, β, x) + logbeta(α, β)
        @test betalogulikelihood(α, β, x) == betalogpdf(α, β, x)

        # Binomial distribution
        n = 9
        p = 0.4
        k = 6
        @test binomlogupdf(n, p, k) ≈ binomlogpdf(n, p, k) + log(n + 1)
        @test binomlogulikelihood(n, p, k) == binomlogpdf(n, p, k)

        # Chi-squared distribution
        k = 4.2
        x = 3.1
        @test chisqlogupdf(k, x) ≈ chisqlogpdf(k, x) + k / 2 * log(2) + loggamma(k / 2)
        @test chisqlogulikelihood(k, x) ≈ chisqlogpdf(k, x) + log(x) + x / 2

        # F distribution
        ν1 = 0.9
        ν2 = 1.5
        x = 2.1
        @test fdistlogupdf(ν1, ν2, x) ≈ fdistlogpdf(ν1, ν2, x) + logbeta(ν1 / 2, ν2 / 2) - (ν1 * log(ν1) + ν2 * log(ν2)) / 2
        @test fdistlogulikelihood(ν1, ν2, x) ≈ fdistlogpdf(ν1, ν2, x) - (ν1 / 2 - 1) * log(x)

        # Gamma distribution
        k = 1.4
        θ = 2.3
        x = 1.9
        @test gammalogupdf(k, θ, x) ≈ gammalogpdf(k, θ, x) + loggamma(k) + k * log(θ)
        @test gammalogulikelihood(k, θ, x) ≈ gammalogpdf(k, θ, x) + log(x)
    end

    @testset "fallback" begin
        # Hyper-geometric distribution
        ms = 2
        mf = 3
        n = 4
        x = 2
        @test hyperlogupdf(ms, mf, n, x) == hyperlogpdf(ms, mf, n, x)
        @test hyperlogulikelihood(ms, mf, n, x) == hyperlogpdf(ms, mf, n, x)

        # Non-central beta distribution
        α = 0.8
        β = 2.1
        λ = 1.1
        x = 0.8
        @test nbetalogupdf(α, β, λ, x) == nbetalogpdf(α, β, λ, x)
        @test nbetalogulikelihood(α, β, λ, x) == nbetalogpdf(α, β, λ, x)

        # Negative binomial distribution
        r = 3
        p = 0.7
        x = 2
        @test nbinomlogupdf(r, p, x) == nbinomlogpdf(r, p, x)
        @test nbinomlogulikelihood(r, p, x) == nbinomlogpdf(r, p, x)
    end
end
