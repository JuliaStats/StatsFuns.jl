tests = ["basicfuns", "rmath", "generic", "binomial"]

for t in tests
    fp = "$t.jl"
    println("* running $fp")
    include(fp)
end
