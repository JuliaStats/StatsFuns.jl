tests = ["rmath", "generic", "misc", "chainrules", "inverse"]

for t in tests
    fp = "$t.jl"
    println("* running $fp")
    include(fp)
end
