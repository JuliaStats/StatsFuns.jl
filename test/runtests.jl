tests = ["basicfuns", "rmath"]

for t in tests
    to_test = true
    # if VERSION < v"0.4.0-dev"
    #     to_test = false
    # end

    if to_test
        fp = "$t.jl"
        println("* running $fp")
        include(fp)
    end
end
