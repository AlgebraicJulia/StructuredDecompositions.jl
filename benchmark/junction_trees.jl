names = ("mycielskian4", "dwt_59", "can_292", "lshp3466", "wing", "144", "333SP")
ssmc = ssmc_db()

for name in names
    # download graph
    graph = mmread(joinpath(fetch_ssmc(ssmc[ssmc.name .== name, :], format="MM")[1], "$(name).mtx"))

    # remove self loops
    fkeep!((i, j, v) -> i != j, graph)

    # construct permutation
    order, index = permutation(graph, AMD())

    # construct benchmarks
    SUITE["junction trees"][name]["StructuredDecompositions"] = @benchmarkable eliminationgraph($graph; alg=$order)
    SUITE["junction trees"][name]["QDLDL"] = @benchmarkable QDLDL.qdldl($(graph + (1.0)I); perm=$order, logical=true)
end
