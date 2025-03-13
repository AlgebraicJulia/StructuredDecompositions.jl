using Documenter
using Literate

@info "Loading StructuredDecompositions"
using StructuredDecompositions

# Set Literate.jl config if not being compiled on recognized service.
config = Dict{String,String}()
if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
  config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/AlgebraicJulia/StructuredDecompositions.jl/blob/gh-pages/dev"
  config["repo_root_url"] = "https://github.com/AlgebraicJulia/StructuredDecompositions.jl/blob/master/docs"
end

# This block of code will compile out examples found in an `examples/` directory in the root and compile them to docs pages
# const literate_dir = joinpath(@__DIR__, "..", "examples")
# const generated_dir = joinpath(@__DIR__, "src", "examples")
#
# for (root, dirs, files) in walkdir(literate_dir)
#   out_dir = joinpath(generated_dir, relpath(root, literate_dir))
#   for file in files
#     f, l = splitext(file)
#     if l == ".jl" && !startswith(f, "_")
#       Literate.markdown(joinpath(root, file), out_dir;
#         config=config, documenter=true, credit=false)
#       Literate.notebook(joinpath(root, file), out_dir;
#         execute=true, documenter=true, credit=false)
#     end
#   end
# end

@info "Building Documenter.jl docs"
makedocs(
  modules = [StructuredDecompositions],
  format = Documenter.HTML(),
  sitename = "StructuredDecompositions.jl",
  doctest = false,
  checkdocs = :none,
  pages = [
    "StructuredDecompositions.jl" => "index.md",
    "Decompositions" => "pages/Decompositions.md",
    "DecidingSheaves" => "pages/DecidingSheaves.md",
    "FunctorUtils" => "pages/FunctorUtils.md",
    "API" => [
        "Decompositions" => "api/Decompositions.md",
        "DecidingSheaves" => "api/DecidingSheaves.md",
        "FunctorUtils" => "api/FunctorUtils.md",
    ]
  ]
)

@info "Deploying docs"
deploydocs(
  target = "build",
  repo = "github.com/AlgebraicJulia/StructuredDecompositions.jl.git",
  branch = "gh-pages"
)
