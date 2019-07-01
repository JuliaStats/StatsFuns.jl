tests = ["basicfuns", "rmath", "generic", "misc"]

for t in tests
    fp = "$t.jl"
    println("* running $fp")
    include(fp)
end
