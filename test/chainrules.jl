using StatsFuns, Test
using ChainRulesCore
using ChainRulesTestUtils
using Random

@testset "chainrules" begin
    x = exp(randn())
    y = exp(randn())
    z = logistic(randn())
    test_frule(betalogpdf, x, y, z)
    test_rrule(betalogpdf, x, y, z)

    x = exp(randn())
    y = exp(randn())
    z = exp(randn())
    test_frule(gammalogpdf, x, y, z)
    test_rrule(gammalogpdf, x, y, z)

    x = exp(randn())
    y = exp(randn())
    test_frule(chisqlogpdf, x, y)
    test_rrule(chisqlogpdf, x, y)

    x = exp(randn())
    y = exp(randn())
    z = exp(randn())
    test_frule(fdistlogpdf, x, y, z)
    test_rrule(fdistlogpdf, x, y, z)

    x = exp(randn())
    y = randn()
    test_frule(tdistlogpdf, x, y)
    test_rrule(tdistlogpdf, x, y)

    x = rand(1:100)
    y = logistic(randn())
    z = rand(1:x)
    test_frule(binomlogpdf, x, y, z)
    test_rrule(binomlogpdf, x, y, z)

    x = exp(randn())
    y = rand(1:100)
    test_frule(poislogpdf, x, y)
    test_rrule(poislogpdf, x, y)

    # test special case λ = 0
    _, pb = rrule(poislogpdf, 0.0, 0)
    _, x̄1, _ = pb(1)
    @test x̄1 == -1
    _, pb = rrule(poislogpdf, 0.0, 1)
    _, x̄1, _ = pb(1)
    @test x̄1 == Inf
end
