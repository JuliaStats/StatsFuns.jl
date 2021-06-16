tests = ["rmath", "generic", "misc","MVN"]

for t in tests
    fp = "$t.jl"
    println("* running $fp")
    include(fp)
end
