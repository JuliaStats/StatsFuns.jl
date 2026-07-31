using StatsFuns
using ForwardDiff: Dual
using Test

# `normpdf` and `normlogpdf` carry `@inline` (#223). Without it they are not
# inlined for `ForwardDiff.Dual` arguments, and a loop calling them costs
# noticeably more than a loop calling `log` on a `Dual`. Measuring against that
# baseline instead of against absolute times keeps the bounds independent of the
# machine the tests run on.
#
# Only these two functions are checked: for the other functions annotated in
# #223 the difference between the inlined and the non-inlined version is a
# factor 1.04-1.33, which is too small to assert on shared CI runners.

function accumulate_f(f, μ, σ, xs)
    s = zero(μ)
    for x in xs
        s += f(μ, σ, x)
    end
    return s
end

function accumulate_log(μ, _, xs)
    s = zero(μ)
    for x in xs
        s += log(μ + x)
    end
    return s
end

function mintime(f, args...)
    t = typemax(UInt64)
    for _ in 1:50
        t0 = time_ns()
        f(args...)
        t = min(t, time_ns() - t0)
    end
    return t
end

@testset "Inlining for Dual arguments" begin
    μ = Dual{:StatsFunsInlining}(0.1, 1.0, 0.0, 0.0, 0.0)
    σ = 1.3
    xs = collect(range(4.0, 9.0; length = 4096))  # positive, so `log` stays in domain

    bounds = ((normpdf, 2.0), (normlogpdf, 1.4))

    # compile everything before timing anything
    accumulate_log(μ, σ, xs)
    for (f, _) in bounds
        accumulate_f(f, μ, σ, xs)
    end

    # interleave baseline and target measurements so that a load spike during one
    # of them cannot skew the ratio in only one direction
    baseline = typemax(UInt64)
    times = Dict{Any, UInt64}(f => typemax(UInt64) for (f, _) in bounds)
    for _ in 1:2
        baseline = min(baseline, mintime(accumulate_log, μ, σ, xs))
        for (f, _) in bounds
            times[f] = min(times[f], mintime(accumulate_f, f, μ, σ, xs))
        end
    end

    @testset "$f" for (f, bound) in bounds
        @test times[f] / baseline < bound
    end
end
