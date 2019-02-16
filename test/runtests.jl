tests = ["basicfuns", "rmath", "generic","tvpack",]

for t in tests
    fp = "$t.jl"
    println("* running $fp")
    include(fp)
end
