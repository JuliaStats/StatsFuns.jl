using StatsFuns, Test
using ChainRulesCore
using ChainRulesTestUtils
using Random

@testset "chainrules" begin
    # `test_frule`/`test_rrule` finite-difference the log-pdfs, which throw `DomainError`
    # outside their domain. ChainRulesTestUtils evaluates the reference at `input ± reach`
    # with `reach = max_range * |tangent|`. With its defaults (`max_range = 1e-2` and, for
    # `Float64`, `|tangent| ≤ 9` since tangents are drawn from `-9:0.01:9`) the reach is at
    # most `0.09`. Drawing every differentiated argument ≥ 0.1 from its boundary therefore
    # keeps the stencil inside the domain by construction (`0.09 < 0.1`), for any draw or
    # seed, so we can rely on the CRTU defaults without a custom `fdm` or tangents. This
    # margin assumes those CRTU defaults; if `max_range` or the tangent range grows, the
    # `0.1` anchor below must grow with it. A fixed seed makes any failure reproducible.
    Random.seed!(1234)
    pos() = 0.1 + randexp()    # positive args (shape/df/rate/λ/eval point)
    unit01() = 0.1 + 0.8rand() # args in [0.1, 0.9): beta `x`, binom `p`

    x = pos()
    y = pos()
    z = unit01()
    test_frule(betalogpdf, x, y, z)
    test_rrule(betalogpdf, x, y, z)

    x = pos()
    y = pos()
    z = pos()
    test_frule(gammalogpdf, x, y, z)
    test_rrule(gammalogpdf, x, y, z)

    x = pos()
    y = pos()
    test_frule(chisqlogpdf, x, y)
    test_rrule(chisqlogpdf, x, y)

    x = pos()
    y = pos()
    z = pos()
    test_frule(fdistlogpdf, x, y, z)
    test_rrule(fdistlogpdf, x, y, z)

    x = pos()
    y = randn()  # location: unbounded
    test_frule(tdistlogpdf, x, y)
    test_rrule(tdistlogpdf, x, y)

    x = rand(1:100)  # `n`: integer, NoTangent
    y = unit01()
    z = rand(1:x)    # `k`: integer, NoTangent
    test_frule(binomlogpdf, x, y, z)
    test_rrule(binomlogpdf, x, y, z)

    x = pos()
    y = rand(1:100)  # `k`: integer, NoTangent
    test_frule(poislogpdf, x, y)
    test_rrule(poislogpdf, x, y)

    # test special case λ = 0
    _, pb = rrule(poislogpdf, 0.0, 0)
    _, x̄1, _ = pb(1)
    @test x̄1 == -1
    _, pb = rrule(poislogpdf, 0.0, 1)
    _, x̄1, _ = pb(1)
    @test x̄1 == Inf

    h = randn()
    a = randn()
    test_frule(owens_t, h, a)
    test_rrule(owens_t, h, a)
end
