tests = ["basicfuns", "rmath", "generic", "misc", "chainrules"]

for t in tests
    fp = "$t.jl"
    println("* running $fp")
    include(fp)
end
