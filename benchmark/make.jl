using BenchmarkTools: BenchmarkGroup, prettytime, prettymemory
using PkgBenchmark: benchmarkpkg
using StructuredDecompositions


const DEFAULT_NAME = "benchmarks.md"


function writemarkdown(io::IO, group::BenchmarkGroup)
    println(io, "# Benchmarks")
    println(io)

    for colors in (2, 3)
        println(io, "## $colors Coloring")
        println(io)
        println(io, "| library | bags | time | gctime | memory | allocs |")
        println(io, "| :------ | :--- | :----| :----- | :----- | :----- |")

        trial = group["graph coloring fixed"]["$colors coloring"]["StructuredDecompositions"]["1 bag"]
        estimate = minimum(trial)
        println(io, "| StructuredDecompositions | 1 | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |")

        for bags in (2, 4, 8, 12)
            trial = group["graph coloring fixed"]["$colors coloring"]["StructuredDecompositions"]["$bags bags"]
            estimate = minimum(trial)
            println(io, "| StructuredDecompositions | $bags | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |")
        end

        trial = group["graph coloring fixed"]["$colors coloring"]["HomSearch"]
        estimate = minimum(trial)
        println(io, "| HomSearch |  | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |")

        trial = group["graph coloring fixed"]["$colors coloring"]["SimpleGraphAlgorithms"]
        estimate = minimum(trial)
        println(io, "| SimpleGraphAlgorithms |  | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |")
        println(io)
    end

    for colors in (4,)
        println(io, "## $colors Coloring")
        println(io)
        println(io, "| library | bags | time | gctime | memory | allocs |")
        println(io, "| :------ | :--- | :--- | :----- | :----- | :----- |")
        
        for bags in (4, 8, 12)
            trial = group["graph coloring fixed"]["$colors coloring"]["StructuredDecompositions"]["$bags bags"]
            estimate = minimum(trial)
            println(io, "| StructuredDecompositions | $bags | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |")
        end

        trial = group["graph coloring fixed"]["4 coloring"]["HomSearch"]
        estimate = minimum(trial)
        println(io, "| HomSearch |  | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |")

        trial = group["graph coloring fixed"]["$colors coloring"]["SimpleGraphAlgorithms"]
        estimate = minimum(trial)
        println(io, "| SimpleGraphAlgorithms |  | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |")
    end
end


function postprocess(group::BenchmarkGroup)
    open(io -> writemarkdown(io, group), get(ARGS, 1, DEFAULT_NAME); write=true)
    group
end


benchmarkpkg(StructuredDecompositions; postprocess)
