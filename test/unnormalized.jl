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

        # Negative binomial distribution
        r = 3
        p = 0.7
        x = 2
        @test nbinomlogupdf(r, p, x) ≈ nbinomlogpdf(r, p, x) - xlogy(r, p)
        @test nbinomlogulikelihood(r, p, x) == nbinomlogpdf(r, p, x)

        # Normal distribution
        x = 3.1
        @test normlogupdf(x) ≈ normlogpdf(x) + log(2 * π) / 2
        @test normlogulikelihood(x) ≈ normlogpdf(x) + log(2 * π) / 2
        for (μ, σ) in ((-0.1, 2.4), (-0.1, 0.0), (x, 0.0))
            @test normlogupdf(μ, σ, x) ≈ (iszero(σ) ? normlogpdf(μ, σ, x) : normlogpdf(μ, σ, x) + log(2 * π) / 2 + log(σ))
            @test normlogulikelihood(μ, σ, x) ≈ normlogpdf(μ, σ, x) + log(2 * π) / 2
        end

        # Poisson distribution
        λ = 3.5
        x = 2
        @test poislogupdf(λ, x) ≈ poislogpdf(λ, x) + λ
        @test poislogulikelihood(λ, x) ≈ poislogpdf(λ, x) + loggamma(x + 1)

        # Student's t distribution
        ν = 4.5
        x = 1.3
        @test tdistlogupdf(ν, x) ≈ tdistlogpdf(ν, x) - loggamma((ν + 1) / 2) + (log(π) + log(ν)) / 2 + loggamma(ν / 2)
        @test tdistlogulikelihood(ν, x) ≈ tdistlogpdf(ν, x) + log(π) / 2

        # Wilcoxon signed rank distribution
        n = 5
        W = 7
        @test signranklogupdf(n, W) ≈ signranklogpdf(n, W) + n * log(2)
        @test signranklogulikelihood(n, W) == signranklogpdf(n, W)

        # Wilcoxon rank sum distribution
        nx = 3
        ny = 4
        U = 5
        @test wilcoxlogupdf(nx, ny, U) ≈ wilcoxlogpdf(nx, ny, U) + first(logabsbinomial(nx + ny, nx))
        @test wilcoxlogulikelihood(nx, ny, U) == wilcoxlogpdf(nx, ny, U)
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

        # Non-central chi-squared distribution
        k = 3.0
        λ = 1.5
        x = 4.2
        @test nchisqlogupdf(k, λ, x) == nchisqlogpdf(k, λ, x)
        @test nchisqlogulikelihood(k, λ, x) == nchisqlogpdf(k, λ, x)

        # Non-central F distribution
        k1 = 2.0
        k2 = 3.0
        λ = 1.0
        x = 1.5
        @test nfdistlogupdf(k1, k2, λ, x) == nfdistlogpdf(k1, k2, λ, x)
        @test nfdistlogulikelihood(k1, k2, λ, x) == nfdistlogpdf(k1, k2, λ, x)

        # Non-central t distribution
        k = 5.0
        λ = 1.2
        x = 2.0
        @test ntdistlogupdf(k, λ, x) == ntdistlogpdf(k, λ, x)
        @test ntdistlogulikelihood(k, λ, x) == ntdistlogpdf(k, λ, x)
    end
end
