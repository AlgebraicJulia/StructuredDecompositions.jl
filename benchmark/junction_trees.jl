names = ("mycielskian4", "can_292", "wing", "333SP")
ssmc = ssmc_db()

for name in names
    # download graph
    graph = mmread(joinpath(fetch_ssmc(ssmc[ssmc.name .== name, :], format="MM")[1], "$(name).mtx"))

    # remove self loops
    fkeep!((i, j, v) -> i != j, graph)

    # construct permutation
    order, index = permutation(graph, AMD())

    # construct benchmarks
    SUITE["junction trees"]["StructuredDecompositions"][name]["maximal"] = @benchmarkable junctiontree($graph; alg=$order, snd=$(Maximal()))
    SUITE["junction trees"]["StructuredDecompositions"][name]["fundamental"] = @benchmarkable junctiontree($graph; alg=$order, snd=$(Fundamental()))
    SUITE["junction trees"]["StructuredDecompositions"][name]["nodal"] = @benchmarkable junctiontree($graph; alg=$order, snd=$(Nodal()))
    SUITE["junction trees"]["QDLDL"][name] = @benchmarkable qdldl($(graph + (1.0)I); perm=$order, logical=true)
end
