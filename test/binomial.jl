using StatsFuns
using Base.Test

N = (1, 1, 8, 20, 20, 20);
P = (0.5, 0.7, 0.6, 0.34, 0.89, 0.53);

println("\ttesting binompdf ...");

for i = 1:length(N)
    for x = Int64(0):N[i]
    bin = factorial(Int64(N[i]))/factorial(x)/factorial(N[i] - x);
    v = bin*(P[i]^x)*(1-P[i])^(N[i] - x);
    @test ≈(v, binompdf(N[i], P[i], x), atol = 1e-15)
    end
end

println("\ttesting binomlogpdf ...");
for i = 1:length(N)
    for x = Int64(0):N[i]
    bin = factorial(Int64(N[i]))/factorial(x)/factorial(N[i] - x);
    v = bin*(P[i]^x)*(1-P[i])^(N[i] - x);
    @test ≈(log(v), binomlogpdf(N[i], P[i], x), atol = 1e-15)
    end
end
