using StatsFuns
using Base.Test

#Test Compares the difference in values w.r.t Rmath implementation

N = Int64.( (0, 1, 1, 8, 20, 20, 20, 150, 700, 900, 1000) );
P = (0.0, 0.5, 0.7, 0.6, 0.34, 0.89, 0.53, 0.77, 0.98, 0.5, 0.29);

@testset "binompdf" begin
for i = 1:length(N)
    for x = Int64(0):N[i]
    if isfinite(binompdf(N[i], P[i], x))
        @test abs(binompdf(N[i], P[i], x) - Rmath.dbinom(x, N[i], P[i])) < 1e-15
    end
end
end
end

@testset "binomlogpdf" begin
for i = 1:length(N)
    for x = Int64(0):N[i]
        val = abs(binomlogpdf(N[i], P[i], x) - Rmath.dbinom(x, N[i], P[i], true)/binomlogpdf(N[i], P[i], x)
        if isfinite(val) && binomlogpdf(N[i], P[i], x) != Rmath.dbinom(x, N[i], P[i], true)
            @test  val < 2.24e-14
        end
    end
end
end
