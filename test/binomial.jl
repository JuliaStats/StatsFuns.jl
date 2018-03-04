using StatsFuns
using Base.Test

#Test Compares the difference in values w.r.t Rmath implementation

N = (1, 1, 8, 20, 20, 20, 150, 700, 900, 1000);
P = (0.5, 0.7, 0.6, 0.34, 0.89, 0.53, 0.77, 0.98, 0.5, 0.29);

println("\ttesting binompdf ...");

for i = 1:length(N)
    for x = Int64(0):N[i]
@test abs(Rmath.dbinom(x, N[i], P[i])/binompdf(N[i], P[i], x) - 1.0) < 1e-14
    end
end

println("\ttesting binomlogpdf ...");
for i = 1:length(N)
    for x = Int64(0):N[i]
    @test  abs(Rmath.dbinom(x, N[i], P[i], true)/binomlogpdf(N[i], P[i], x) - 1.0) < 1e-14
    end
end
