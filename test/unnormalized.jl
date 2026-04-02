using StatsFuns
using SpecialFunctions
using Test

@testset "logupdf + logulikelihood" begin
    @testset "optimized" begin
        # Beta distribution
        for ((α, β, x), T) in (((0.1, 2.5, 0.7), Float64), ((1, 2.5, 0.7), Float64), ((0.1f0, 2.5f0, 0.7f0), Float32), ((1, 2.5f0, 0.7f0), Float32))
            @test (@inferred betalogupdf(α, β, x))::T ≈ betalogpdf(α, β, x) + logbeta(α, β)
            @test (@inferred betalogulikelihood(α, β, x))::T == betalogpdf(α, β, x)
        end

        # Binomial distribution
        for ((n, p, k), T) in (((9, 0.4, 6), Float64), ((9.0f0, 0.4f0, 6.0f0), Float32), ((9, 0.4f0, 6), Float32))
            @test (@inferred binomlogupdf(n, p, k))::T ≈ binomlogpdf(n, p, k) + log(T(n) + 1)
            @test (@inferred binomlogulikelihood(n, p, k))::T == binomlogpdf(n, p, k)
        end

        # Chi-squared distribution
        for ((k, x), T) in (((4.2, 3.1), Float64), ((4, 3.1), Float64), ((4.2f0, 3.1f0), Float32), ((4, 3.1f0), Float32))
            @test (@inferred chisqlogupdf(k, x))::T ≈ chisqlogpdf(k, x) + T(k) / 2 * log(T(2)) + loggamma(T(k) / 2)
            @test (@inferred chisqlogulikelihood(k, x))::T ≈ chisqlogpdf(k, x) + log(T(x)) + T(x) / 2
        end

        # F distribution
        for ((ν1, ν2, x), T) in (((0.9, 1.5, 2.1), Float64), ((1, 2, 2.1), Float64), ((0.9f0, 1.5f0, 2.1f0), Float32), ((1, 2, 2.1f0), Float32))
            @test (@inferred fdistlogupdf(ν1, ν2, x))::T ≈ fdistlogpdf(ν1, ν2, x) + logbeta(T(ν1) / 2, T(ν2) / 2) - (T(ν1) * log(T(ν1)) + T(ν2) * log(T(ν2))) / 2
            @test (@inferred fdistlogulikelihood(ν1, ν2, x))::T ≈ fdistlogpdf(ν1, ν2, x) - (T(ν1) / 2 - 1) * log(T(x))
        end

        # Gamma distribution
        for ((k, θ, x), T) in (((1.4, 2.3, 1.9), Float64), ((2, 2.3, 1.9), Float64), ((1.4f0, 2.3f0, 1.9f0), Float32), ((2, 2.3f0, 1.9f0), Float32))
            @test (@inferred gammalogupdf(k, θ, x))::T ≈ gammalogpdf(k, θ, x) + loggamma(T(k)) + T(k) * log(T(θ))
            @test (@inferred gammalogulikelihood(k, θ, x))::T ≈ gammalogpdf(k, θ, x) + log(T(x))
        end

        # Negative binomial distribution
        for ((r, p, x), T) in (((3, 0.7, 2), Float64), ((3.0f0, 0.7f0, 2.0f0), Float32), ((3, 0.7f0, 2), Float32))
            @test (@inferred nbinomlogupdf(r, p, x))::T ≈ nbinomlogpdf(r, p, x) - xlogy(T(r), T(p))
            @test (@inferred nbinomlogulikelihood(r, p, x))::T == nbinomlogpdf(r, p, x)
        end
        for ((r, p, x), T) in (((3, 1.5, 2), Float64), ((-1, 0.5, 2), Float64), ((3.0f0, 1.5f0, 2.0f0), Float32))
            @test isnan((@inferred nbinomlogupdf(r, p, x))::T)
        end

        # Normal distribution
        for (z, T) in ((3.1, Float64), (3, Float64), (3.1f0, Float32))
            @test (@inferred normlogupdf(z))::T ≈ normlogpdf(z) + log(T(2π)) / 2
            @test (@inferred normlogulikelihood(z))::T ≈ normlogpdf(z) + log(T(2π)) / 2
        end
        for ((μ, σ, z), T) in (
                ((-0.1, 2.4, 3.1), Float64), ((-0.1, 0.0, 3.1), Float64), ((3.1, 0.0, 3.1), Float64),
                ((-0.1f0, 2.4f0, 3.1f0), Float32), ((-0.1f0, 0.0f0, 3.1f0), Float32), ((3.1f0, 0.0f0, 3.1f0), Float32), ((0, 2.4f0, 3.1f0), Float32),
            )
            @test (@inferred normlogupdf(μ, σ, z))::T ≈ (iszero(σ) ? normlogpdf(μ, σ, z) : normlogpdf(μ, σ, z) + log(T(2π)) / 2 + log(σ))
            @test (@inferred normlogulikelihood(μ, σ, z))::T ≈ normlogpdf(μ, σ, z) + log(T(2π)) / 2
        end

        # Poisson distribution
        for ((λ, x), T) in (((3.5, 2), Float64), ((3.5f0, 2.0f0), Float32), ((3.5f0, 2), Float32))
            @test (@inferred poislogupdf(λ, x))::T ≈ poislogpdf(λ, x) + T(λ)
            @test (@inferred poislogulikelihood(λ, x))::T ≈ poislogpdf(λ, x) + loggamma(T(x) + 1)
        end

        # Student's t distribution
        for ((ν, x), T) in (((4.5, 1.3), Float64), ((5, 1.3), Float64), ((4.5f0, 1.3f0), Float32), ((5, 1.3f0), Float32))
            @test (@inferred tdistlogupdf(ν, x))::T ≈ tdistlogpdf(ν, x) - loggamma((T(ν) + 1) / 2) + (log(T(π)) + log(T(ν))) / 2 + loggamma(T(ν) / 2)
            @test (@inferred tdistlogulikelihood(ν, x))::T ≈ tdistlogpdf(ν, x) + log(T(π)) / 2
        end

        # Wilcoxon signed rank distribution
        for (n, W) in ((5, 7), (5, 7.0))
            @test (@inferred signranklogupdf(n, W))::Float64 ≈ signranklogpdf(n, W) + n * log(2)
            @test (@inferred signranklogulikelihood(n, W))::Float64 == signranklogpdf(n, W)
        end

        # Wilcoxon rank sum distribution
        for (nx, ny, U) in ((3, 4, 5), (3, 4, 5.0))
            @test (@inferred wilcoxlogupdf(nx, ny, U))::Float64 ≈ wilcoxlogpdf(nx, ny, U) + first(logabsbinomial(nx + ny, nx))
            @test (@inferred wilcoxlogulikelihood(nx, ny, U))::Float64 == wilcoxlogpdf(nx, ny, U)
        end
    end

    @testset "fallback" begin
        # Hyper-geometric distribution
        for (ms, mf, n, x) in ((2, 3, 4, 2), (2, 3, 4, 2.0))
            @test (@inferred hyperlogupdf(ms, mf, n, x))::Float64 == hyperlogpdf(ms, mf, n, x)
            @test (@inferred hyperlogulikelihood(ms, mf, n, x))::Float64 == hyperlogpdf(ms, mf, n, x)
        end

        # Non-central beta distribution
        for ((α, β, λ, x), T) in (((0.8, 2.1, 1.1, 0.8), Float64), ((1, 2.1, 1.1, 0.8), Float64), ((0.8f0, 2.1f0, 1.1f0, 0.8f0), Float32), ((1, 2.1f0, 1.1f0, 0.8f0), Float32))
            @test (@inferred nbetalogupdf(α, β, λ, x))::T == nbetalogpdf(α, β, λ, x)
            @test (@inferred nbetalogulikelihood(α, β, λ, x))::T == nbetalogpdf(α, β, λ, x)
        end

        # Non-central chi-squared distribution
        for ((k, λ, x), T) in (((3.0, 1.5, 4.2), Float64), ((3, 1.5, 4.2), Float64), ((3.0f0, 1.5f0, 4.2f0), Float32), ((3, 1.5f0, 4.2f0), Float32))
            @test (@inferred nchisqlogupdf(k, λ, x))::T == nchisqlogpdf(k, λ, x)
            @test (@inferred nchisqlogulikelihood(k, λ, x))::T == nchisqlogpdf(k, λ, x)
        end

        # Non-central F distribution
        for ((k1, k2, λ, x), T) in (((2.0, 3.0, 1.0, 1.5), Float64), ((2, 3, 1.0, 1.5), Float64), ((2.0f0, 3.0f0, 1.0f0, 1.5f0), Float32), ((2, 3, 1.0f0, 1.5f0), Float32))
            @test (@inferred nfdistlogupdf(k1, k2, λ, x))::T == nfdistlogpdf(k1, k2, λ, x)
            @test (@inferred nfdistlogulikelihood(k1, k2, λ, x))::T == nfdistlogpdf(k1, k2, λ, x)
        end

        # Non-central t distribution
        for ((k, λ, x), T) in (((5.0, 1.2, 2.0), Float64), ((5, 1.2, 2.0), Float64), ((5.0f0, 1.2f0, 2.0f0), Float32), ((5, 1.2f0, 2.0f0), Float32))
            @test (@inferred ntdistlogupdf(k, λ, x))::T == ntdistlogpdf(k, λ, x)
            @test (@inferred ntdistlogulikelihood(k, λ, x))::T == ntdistlogpdf(k, λ, x)
        end
    end
end
