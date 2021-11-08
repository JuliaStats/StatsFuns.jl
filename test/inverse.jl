using StatsFuns
using Test
using InverseFunctions

@testset "inverse" begin
    x = 0.7
    InverseFunctions.test_inverse(norminvcdf, x)
    InverseFunctions.test_inverse(norminvccdf, x)
    InverseFunctions.test_inverse(norminvlogcdf, log(x))
    InverseFunctions.test_inverse(norminvlogccdf, log(x))
end
