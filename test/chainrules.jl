using StatsFuns, Test
using ChainRulesTestUtils
using Random

# move upstream?
ChainRulesTestUtils.rand_tangent(rng::Random.AbstractRNG, x::BigFloat) = big(randn(rng))

@testset "chainrules" begin
    x, Δx, x̄ = randn(3)
    y, Δy, ȳ = randn(3)
    z, Δz, z̄ = randn(3)
    Δu = randn()

    x̃ = exp(x)
    ỹ = exp(y)
    z̃ = logistic(z)
    frule_test(betalogpdf, (x̃, Δx), (ỹ, Δy), (z̃, Δz))
    rrule_test(betalogpdf, Δu, (x̃, x̄), (ỹ, ȳ), (z̃, z̄))

    x̃ = exp(x)
    ỹ = exp(y)
    z̃ = exp(z)
    frule_test(gammalogpdf, (x̃, Δx), (ỹ, Δy), (z̃, Δz))
    rrule_test(gammalogpdf, Δu, (x̃, x̄), (ỹ, ȳ), (z̃, z̄))

    x̃ = exp(x)
    ỹ = exp(y)
    z̃ = exp(z)
    frule_test(chisqlogpdf, (x̃, Δx), (ỹ, Δy))
    rrule_test(chisqlogpdf, Δu, (x̃, x̄), (ỹ, ȳ))

    x̃ = exp(x)
    ỹ = exp(y)
    z̃ = exp(z)
    frule_test(fdistlogpdf, (x̃, Δx), (ỹ, Δy), (z̃, Δz))
    rrule_test(fdistlogpdf, Δu, (x̃, x̄), (ỹ, ȳ), (z̃, z̄))

    x̃ = exp(x)
    frule_test(tdistlogpdf, (x̃, Δx), (y, Δy))
    rrule_test(tdistlogpdf, Δu, (x̃, x̄), (y, ȳ))

    # use `BigFloat` to avoid Rmath implementation in finite differencing check
    # (returns `NaN` for non-integer values)
    x̃ = BigFloat(rand(1:100))
    ỹ = big(logistic(y))
    z̃ = BigFloat(rand(1:x̃))
    frule_test(binomlogpdf, (x̃, big(Δx)), (ỹ, big(Δy)), (z̃, big(Δz)))
    rrule_test(binomlogpdf, big(Δu), (x̃, big(x̄)), (ỹ, big(ȳ)), (z̃, big(z̄)))

    x̃ = big(exp(x))
    ỹ = BigFloat(rand(1:100))
    frule_test(poislogpdf, (x̃, big(Δx)), (ỹ, big(Δy)))
    rrule_test(poislogpdf, big(Δu), (x̃, big(x̄)), (ỹ, big(ȳ)))
end
