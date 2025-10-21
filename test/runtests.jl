tests = ["rmath", "generic", "misc", "chainrules", "inverse", "tvpack", "qa"]

for t in tests
    fp = "$t.jl"
    println("* running $fp")
    include(fp)
end
