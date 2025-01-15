using BenchmarkTools: BenchmarkGroup, prettytime, prettymemory
using PkgBenchmark: benchmarkpkg
using StructuredDecompositions


const DEFAULT_NAME = "benchmarks.md"
const COLORS = (2, 3)
const BAGS = (1, 2, 4, 8, 12)


const GRAPHS = (
    (name="mycielskian2", nv=2,       ne=1),
    (name="mycielskian4", nv=11,      ne=23),
    (name="dwt_59",       nv=59,      ne=104),
    (name="can_292",      nv=292,     ne=1124),
    (name="lshp3466",     nv=3466,    ne=10215),
    (name="wing",         nv=62032,   ne=121544),
    (name="144",          nv=144649,  ne=1074393),
    (name="333SP",        nv=3712815, ne=11108633),)


function writemd(io::IO, group::BenchmarkGroup)
    println(io, "# Benchmarks")
    println(io)
    println(io, "To regenerate this file, navigate to the ``benchmarks`` directory and run the command ``julia --project make.jl``.")
    println(io)
    println(io, "## Junction Tree Construction")
    println(io)
    println(io, "| library | name | vertices | edges | time | gctime | memory | allocs |")
    println(io, "| :------ | :--- | :--------| :-----| :----| :----- | :----- | :----- |")

    for graph in GRAPHS
        name = graph[:name]
        nv = graph[:nv]
        ne = graph[:ne]

        for library in ("StructuredDecompositions", "QDLDL")
            trial = group["junction trees"][library][name]
            estimate = minimum(trial)
            println(io, "| $library | $name | $nv | $ne | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |")
        end
    end

    println(io)
    println(io, "## Vertex Coloring")
    println(io)
    println(io, "| library | colors | bags | time | gctime | memory | allocs |")
    println(io, "| :------ | :----- | :--- | :----| :----- | :----- | :----- |")
   
    for nc in COLORS
        library = "StructuredDecompositions"

        for nb in BAGS
            trial = group["graph coloring fixed"]["$nc coloring"][library]["$nb bags"]
            estimate = minimum(trial)
            println(io, "| $library | $nc | $nb | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |")
        end

        for library in ("Catlab", "SimpleGraphAlgorithms")
            trial = group["graph coloring fixed"]["$nc coloring"][library]
            estimate = minimum(trial)
            println(io, "| $library | $nc |     | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |")
        end
    end
end


function postprocess(group::BenchmarkGroup)
    open(io -> writemd(io, group), get(ARGS, 1, DEFAULT_NAME); write=true)
    group
end


benchmarkpkg(StructuredDecompositions; postprocess)
