@testitem "Inverse" begin
    using StatsFuns
    using Test
    using InverseFunctions

    @testset "inverse" begin
        x = 0.7
        InverseFunctions.test_inverse(norminvcdf, x)
        InverseFunctions.test_inverse(norminvccdf, x)
        InverseFunctions.test_inverse(norminvlogcdf, log(x))
        InverseFunctions.test_inverse(norminvlogccdf, log(x))

        InverseFunctions.test_inverse(laplaceinvcdf, x)
        InverseFunctions.test_inverse(laplaceinvccdf, x)
        InverseFunctions.test_inverse(laplaceinvlogcdf, log(x))
        InverseFunctions.test_inverse(laplaceinvlogccdf, log(x))
    end
end
