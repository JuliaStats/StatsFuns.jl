tests = ["basicfuns", "rmath", "generic"]

for t in tests
    fp = "$t.jl"
    println("* running $fp")
    include(fp)
end
