tests = ["owens_t", "rmath", "generic", "misc", "chainrules", "inverse", "inlining", "tvpack", "qa"]

for t in tests
    fp = "$t.jl"
    println("* running $fp")
    include(fp)
end
