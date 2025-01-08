ssmc = ssmc_db()

names = ("598a", "144", "m14b", "auto")
paths = fetch_ssmc(ssmc[in.(ssmc.name, (names,)), :], format="MM")

for name in names
    path = only(fetch_ssmc(ssmc[ssmc.name .== name, :], format="MM"))
    matrix = mmread(joinpath(path, "$(name).mtx"))
    SUITE["JunctionTrees"]["Nodal"][name] = @benchmarkable junctiontree($matrix; snd=Nodal())
    SUITE["JunctionTrees"]["Fundamental"][name] = @benchmarkable junctiontree($matrix; snd=Fundamental())
    SUITE["JunctionTrees"]["Maximal"][name] = @benchmarkable junctiontree($matrix; snd=Maximal())
end
