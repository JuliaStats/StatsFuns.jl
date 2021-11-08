using StatsFuns
using Test
using InverseFunctions

@testset "inverse" begin
    for d in (:norm,)
        @eval begin
            x = 0.7
            InverseFunctions.test_inverse($(Symbol(d,:invcdf)), x)
            InverseFunctions.test_inverse($(Symbol(d,:invccdf)), x)
            InverseFunctions.test_inverse($(Symbol(d,:invlogcdf)), log(x))
            InverseFunctions.test_inverse($(Symbol(d,:invlogccdf)), log(x))
        end
    end
end
