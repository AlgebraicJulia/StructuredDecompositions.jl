using BenchmarkTools: BenchmarkGroup, Trial, prettytime, prettymemory
using Dates
using PkgBenchmark: benchmarkpkg
using StructuredDecompositions


const DEFAULT_NAME = "README.md"
const SUPERNODES = ("maximal", "fundamental", "nodal")
const COLORS = (2, 3)
const BAGS = (1, 2, 4, 8, 12)


const GRAPHS = (
    (name="mycielskian4", nv=11,      ne=23),
    (name="can_292",      nv=292,     ne=1124),
    (name="wing",         nv=62032,   ne=121544),
    (name="333SP",        nv=3712815, ne=11108633),)


function row(trial::Trial)
    estimate = minimum(trial)
    "$(prettytime(estimate.time)) | $(prettymemory(estimate.memory))"
end


function writemd(io::IO, group::BenchmarkGroup)
    println(io, "# Benchmarks")
    println(io)
    println(io, "This file was automatically generated on $(today()). To regenerate it, navigate to the ``benchmark`` directory and run the following command.")
    println(io, "```")
    println(io, "julia --project make.jl")
    println(io, "```")
    println(io)
    println(io, "## Junction Tree Construction")
    println(io)
    println(io, "| library | supernode | name | edges | time | memory |")
    println(io, "| :------ | :-------- | :----| :---- | :--- | :----- |")

    for graph in GRAPHS
        name = graph[:name]
        ne = graph[:ne]
        library = "StructuredDecompositions"

        for snd in SUPERNODES
            trial = group["junction trees"][library][name][snd]
            println(io, "| $library | $snd | $name | $ne | $(row(trial)) |")
        end

        for library in ("QDLDL",)
            trial = group["junction trees"][library][name]
            println(io, "| $library |      | $name | $ne | $(row(trial)) |")
        end
    end

    println(io)
    println(io, "## Vertex Coloring")
    println(io)
    println(io, "| library | bags | colors | time | memory |")
    println(io, "| :------ | :--- | :----- | :----| :----- |")
   
    for nc in COLORS
        library = "StructuredDecompositions"

        for nb in BAGS
            trial = group["graph coloring fixed"]["$nc coloring"][library]["$nb bags"]
            println(io, "| $library | $nb | $nc | $(row(trial)) |")
        end

        for library in ("Catlab", "SimpleGraphAlgorithms")
            trial = group["graph coloring fixed"]["$nc coloring"][library]
            println(io, "| $library |     | $nc | $(row(trial)) |")
        end
    end
end


function postprocess(group::BenchmarkGroup)
    open(io -> writemd(io, group), get(ARGS, 1, DEFAULT_NAME); write=true)
    group
end


benchmarkpkg(StructuredDecompositions; postprocess)
