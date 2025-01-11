using BenchmarkTools: prettytime, prettymemory
using PkgBenchmark
using PkgBenchmark: benchmarkgroup
using StructuredDecompositions


function makebenchmarks(name)
    group = benchmarkgroup(benchmarkpkg(StructuredDecompositions))
    markdown = "# Benchmarks\n\n"

    for colors in (2, 3)
        markdown *= "## $colors Coloring\n\n"
        markdown *= "| library | bags | time | gctime | memory | allocs |\n"
        markdown *= "| :------ | :--- | :----| :----- | :----- | :----- |\n"

        trial = group["graph coloring fixed"]["$colors coloring"]["StructuredDecompositions"]["1 bag"]
        estimate = minimum(trial)
        markdown *= "| StructuredDecompositions | 1 | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |\n"

        for bags in (2, 4, 8, 12)
            trial = group["graph coloring fixed"]["$colors coloring"]["StructuredDecompositions"]["$bags bags"]
            estimate = minimum(trial)
            markdown *= "| StructuredDecompositions | $bags | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |\n"
        end

        trial = group["graph coloring fixed"]["$colors coloring"]["HomSearch"]
        estimate = minimum(trial)
        markdown *= "| HomSearch |  | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |\n\n"
    end


    for colors in (4,)
        markdown *= "## $colors Coloring\n\n"
        markdown *= "| library | bags | time | gctime | memory | allocs |\n"
        markdown *= "| :------ | :--- | :--- | :----- | :----- | :----- |\n"
        
        for bags in (4, 8, 12)
            trial = group["graph coloring fixed"]["$colors coloring"]["StructuredDecompositions"]["$bags bags"]
            estimate = minimum(trial)
            markdown *= "| StructuredDecompositions | $bags | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |\n"
        end

        trial = group["graph coloring fixed"]["4 coloring"]["HomSearch"]
        estimate = minimum(trial)
        markdown *= "| HomSearch |  | $(prettytime(estimate.time)) | $(prettytime(estimate.gctime)) | $(prettymemory(estimate.memory)) | $(estimate.allocs) |\n"
    end

    write(name, markdown)
end


name = isempty(ARGS) ? "bencharks.md" : first(ARGS)
makebenchmarks(name)
