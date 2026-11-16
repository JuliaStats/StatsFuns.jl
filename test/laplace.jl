@testmodule LaplaceRef begin
    using StatsFuns

    # Closed forms in `BigFloat`; Rmath does not provide the Laplace distribution
    logcdf(z) = z < 0 ? z - log(BigFloat(2)) : log1p(-exp(-z) / 2)
    invlogcdf(lp) = lp <= -log(BigFloat(2)) ? lp + log(BigFloat(2)) : -log(-2 * expm1(lp))
    const REFS = (
        laplacepdf => z -> exp(-abs(z)) / 2,
        laplacelogpdf => z -> -abs(z) - log(BigFloat(2)),
        laplacecdf => z -> exp(logcdf(z)),
        laplaceccdf => z -> exp(logcdf(-z)),
        laplacelogcdf => logcdf,
        laplacelogccdf => z -> logcdf(-z),
        laplaceinvcdf => p -> invlogcdf(log(p)),
        laplaceinvccdf => p -> -invlogcdf(log(p)),
        laplaceinvlogcdf => invlogcdf,
        laplaceinvlogccdf => lp -> -invlogcdf(lp),
    )

    # Evaluate the reference with `4p` bits for an argument of precision `p`
    ref(f, x) = setprecision(() -> Dict(REFS)[f](BigFloat(x)), BigFloat, 4 * precision(x))
end

@testitem "Laplace" setup = [LaplaceRef] begin
    using StatsFuns, Test

    @testset "$f(::$T)" for T in (Float16, Float32, Float64, BigFloat), (f, _) in LaplaceRef.REFS
        l2 = T(-log(big(2)))
        c = T(-log(big(2)) / 2)
        xs = if f in (laplaceinvcdf, laplaceinvccdf)
            T[floatmin(T), 1.0e-3, 0.1, 0.25, prevfloat(T(0.5)), 0.5, nextfloat(T(0.5)), 0.75, 0.9, prevfloat(one(T))]
        elseif f in (laplaceinvlogcdf, laplaceinvlogccdf)
            # around the median, where the sign of the result is decided, and around the branch cutoff
            T[
                -1.0e4, -30, -1, nextfloat(l2, -30), prevfloat(l2), l2, nextfloat(l2), nextfloat(l2, 30),
                prevfloat(c), c, nextfloat(c), -0.1, -1.0e-3, -eps(T), -floatmin(T),
            ]
        else
            T[-50, -2, -1, -0.5, -1.0e-3, -floatmin(T), 0, floatmin(T), 1.0e-3, 0.5, 1, 2, 50]
        end
        @testset "x = $x" for x in xs
            y = @inferred f(x)
            @test y isa T
            r = LaplaceRef.ref(f, x)
            @test abs(y - r) <= 3 * eps(T(r))
        end
    end

    @testset "special values" begin
        @test laplacelogcdf(-Inf) === -Inf
        @test laplacelogccdf(Inf) === -Inf
        @test laplaceinvcdf(0.0) === -Inf
        @test laplaceinvcdf(1.0) === Inf
        @test laplaceinvlogcdf(-Inf) === -Inf
        @test laplaceinvlogcdf(0.0) === Inf
        @test laplaceinvlogccdf(0.0) === -Inf
        @testset "$f" for (f, _) in LaplaceRef.REFS
            @test isnan(f(NaN))
        end
        @test_throws DomainError laplaceinvcdf(1.5)
        @test_throws DomainError laplaceinvlogcdf(0.1)
        @test laplaceinvlogcdf(-1) === laplaceinvlogcdf(-1.0)
    end

    @testset "low part of log(2): T = $T" for T in (Float16, Float32, Float64, BigFloat)
        # `x` at the default precision of `T`, outside of `setprecision`
        x = T(logtwo)
        lo = setprecision(() -> big(logtwo) - x, BigFloat, 4 * precision(T))
        @test StatsFuns._logtwo_lo(T)::T == T(lo)
    end

    @testset "location and scale" begin
        μ, θ = 1.5, 2.0
        @testset "x = $x" for x in (-3.0, 1.5, 4.0)
            z = (x - μ) / θ
            @test laplacepdf(μ, θ, x) ≈ LaplaceRef.ref(laplacepdf, z) / θ rtol = 3 * eps()
            @test laplacelogpdf(μ, θ, x) ≈ LaplaceRef.ref(laplacelogpdf, z) - log(θ) rtol = 3 * eps()
            @testset "$f" for f in (laplacecdf, laplaceccdf, laplacelogcdf, laplacelogccdf)
                @test f(μ, θ, x) ≈ LaplaceRef.ref(f, z) rtol = 3 * eps()
            end
        end
        @testset "p = $p" for p in (0.1, 0.5, 0.9)
            @test laplaceinvcdf(μ, θ, p) ≈ μ + θ * LaplaceRef.ref(laplaceinvcdf, p) rtol = 3 * eps()
            @test laplaceinvccdf(μ, θ, p) ≈ μ + θ * LaplaceRef.ref(laplaceinvccdf, p) rtol = 3 * eps()
            @test laplaceinvlogcdf(μ, θ, log(p)) ≈ μ + θ * LaplaceRef.ref(laplaceinvlogcdf, log(p)) rtol = 3 * eps()
            @test laplaceinvlogccdf(μ, θ, log(p)) ≈ μ + θ * LaplaceRef.ref(laplaceinvlogccdf, log(p)) rtol = 3 * eps()
        end

        # parameters are promoted before computing in low precision
        @test laplacelogpdf(0.0, 3.0f0, 1.0) === laplacelogpdf(0.0, 3.0, 1.0)
        @testset "$f" for f in (laplaceinvcdf, laplaceinvccdf)
            @test f(0.0, 1.0, 0.9f0) === f(Float64(0.9f0))
        end
        @testset "$f" for f in (laplaceinvlogcdf, laplaceinvlogccdf)
            @test f(0.0, 1.0, -0.1f0) === f(Float64(-0.1f0))
        end
    end
end

@testitem "tdistinvlogcdf via laplaceinvlogcdf" begin
    using StatsFuns, Test

    # closed forms of the t quantile for ν = 1 (Cauchy) and ν = 2
    function tquantile(ν, logp)
        return setprecision(BigFloat, 1024) do
            p = exp(big(logp))
            ν == 1 ? -cot(big(π) * p) : (2p - 1) / sqrt(2p * (1 - p))
        end
    end
    @testset "T = $T, ν = $ν" for T in (Float32, Float64), ν in (1, 2)
        l2 = T(-log(big(2)))
        @testset "logp = $logp" for logp in T[-3, -1, nextfloat(l2, -3), l2, nextfloat(l2, 3), -0.5, -0.1]
            t = @inferred(tdistinvlogcdf(T(ν), logp))::T
            @test t ≈ tquantile(ν, logp) rtol = 64 * eps(T)
        end
    end
end

@testitem "Laplace ForwardDiff" tags = [:autodiff] begin
    using StatsFuns, Test
    using ForwardDiff: derivative

    # `1` on the lower half, `u / (2 - u)` with `u = exp(-z)` and `p / (1 - p)` otherwise
    @testset "z = $z" for z in (-1.0e10, -1.0, 0.0, 0.5, 3.0, 1.0e2)
        u = exp(-big(z))
        r = z < 0 ? big(1) : u / (2 - u)
        @test derivative(laplacelogcdf, z) ≈ r rtol = 3 * eps()
    end
    @testset "lp = $lp" for lp in (-1.0e10, -1.0, -log(2), -0.5, -1.0e-3, -1.0e-10)
        p = exp(big(lp))
        r = p <= 0.5 ? big(1) : p / -expm1(big(lp))
        @test derivative(laplaceinvlogcdf, lp) ≈ r rtol = 3 * eps()
    end
end
