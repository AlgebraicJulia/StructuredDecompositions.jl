var documenterSearchIndex = {"docs":
[{"location":"api/FunctorUtils/#FunctorUtils","page":"FunctorUtils","title":"FunctorUtils","text":"","category":"section"},{"location":"api/FunctorUtils/","page":"FunctorUtils","title":"FunctorUtils","text":"Modules = [StructuredDecompositions.FunctorUtils]","category":"page"},{"location":"pages/FunctorUtils/#FunctorUtils","page":"FunctorUtils","title":"Functor Utils","text":"","category":"section"},{"location":"pages/FunctorUtils/","page":"FunctorUtils","title":"FunctorUtils","text":"Functor Utils only includes 4 functions and builds off of Decompositions.","category":"page"},{"location":"pages/FunctorUtils/","page":"FunctorUtils","title":"FunctorUtils","text":"We first define the forgetful functors vs.","category":"page"},{"location":"pages/FunctorUtils/","page":"FunctorUtils","title":"FunctorUtils","text":"function vs(X::Graph) FinSet(length(vertices(X))) end\nfunction vs(f::ACSetTransformation) components(f)[1] end","category":"page"},{"location":"pages/FunctorUtils/","page":"FunctorUtils","title":"FunctorUtils","text":"We also define the functor skeleton taking set to the skeleton of the set.","category":"page"},{"location":"pages/FunctorUtils/","page":"FunctorUtils","title":"FunctorUtils","text":"function skeleton(s::FinSet) FinSet(length(s)) end\nfunction skeleton(f::FinFunction)\n  (dd, cc) = (dom(f), codom(f))\n  ℓ = isempty(dd) ? Int[] : [findfirst(item -> item == f(x), collect(cc)) for x ∈ collect(dd)]\n  FinFunction(ℓ, skeleton(dd), skeleton(cc))\nend","category":"page"},{"location":"api/JunctionTrees/#JunctionTrees","page":"JunctionTrees","title":"JunctionTrees","text":"","category":"section"},{"location":"api/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"CurrentModule = StructuredDecompositions.JunctionTrees","category":"page"},{"location":"api/JunctionTrees/#Trees","page":"JunctionTrees","title":"Trees","text":"","category":"section"},{"location":"api/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"Tree\nSinglyLinkedList\neliminationtree\neliminationtree!","category":"page"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.Tree","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.Tree","text":"Tree <: AbstractUnitRange{Int}\n\nA rooted tree. This type implements the indexed tree interface.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.SinglyLinkedList","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.SinglyLinkedList","text":"SinglyLinkedList{Init <: AbstractScalar{Int}, Next <: AbstractVector{Int}}\n\nA singly linked list. This type implements the iteration interface.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.eliminationtree","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.eliminationtree","text":"eliminationtree(matrix::AbstractMatrix;\n    alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)\n\nA non-mutating version of eliminationtree!.\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.eliminationtree!","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.eliminationtree!","text":"eliminationtree!(matrix::AbstractMatrix;\n    alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)\n\nConstruct a tree-depth decomposition of a connected simple graph. See junctiontree! for the meaning of alg.\n\njulia> using StructuredDecompositions\n\njulia> graph = [\n           0 1 1 0 0 0 0 0\n           1 0 1 0 0 1 0 0\n           1 1 0 1 1 0 0 0\n           0 0 1 0 1 0 0 0\n           0 0 1 1 0 0 1 1\n           0 1 0 0 0 0 1 0\n           0 0 0 0 1 1 0 1\n           0 0 0 0 1 0 1 0\n       ];\n\njulia> label, tree = eliminationtree(graph);\n\njulia> tree\n8-element Tree:\n8\n└─ 7\n   ├─ 5\n   │  ├─ 1\n   │  └─ 4\n   │     └─ 3\n   │        └─ 2\n   └─ 6\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#Junction-Trees","page":"JunctionTrees","title":"Junction Trees","text":"","category":"section"},{"location":"api/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"Bag\nJunctionTree\njunctiontree!\njunctiontree\ntreewidth!\ntreewidth\nseparator\nresidual\nrelative","category":"page"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.Bag","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.Bag","text":"Bag <: AbstractVector{Int}\n\nA bag of a junction tree.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.JunctionTree","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.JunctionTree","text":"JunctionTree <: AbstractVector{Bag}\n\nA junction tree. This type implements the indexed tree interface.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.junctiontree!","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.junctiontree!","text":"junctiontree!(matrix::SparseMatrixCSC;\n    alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,\n    snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)\n\nConstruct a tree decomposition of a connected simple graph. The vertices of the graph are first ordered by a fill-reducing permutation computed by the algorithm alg. The size of the resulting decomposition is determined by the supernode partition snd.\n\njulia> using StructuredDecompositions\n\njulia> graph = [\n           0 1 1 0 0 0 0 0\n           1 0 1 0 0 1 0 0\n           1 1 0 1 1 0 0 0\n           0 0 1 0 1 0 0 0\n           0 0 1 1 0 0 1 1\n           0 1 0 0 0 0 1 0\n           0 0 0 0 1 1 0 1\n           0 0 0 0 1 0 1 0\n       ];\n\njulia> label, tree = junctiontree(graph);\n\njulia> tree\n6-element JunctionTree:\n[6, 7, 8]\n├─ [1, 6, 7]\n├─ [4, 6, 8]\n│  └─ [3, 4, 6]\n│     └─ [2, 3, 6]\n└─ [5, 7, 8]\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.junctiontree","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.junctiontree","text":"junctiontree(matrix::AbstractMatrix;\n    alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,\n    snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)\n\nA non-mutating version of junctiontree!.\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.treewidth!","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.treewidth!","text":"treewidth!(matrix::SparseMatrixCSC;\n    alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)\n\nCompute an upper bound to the tree width of a simple graph. See junctiontree! for the meaning of alg.\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.treewidth","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.treewidth","text":"treewidth(tree::JunctionTree)\n\nCompute the width of a junction tree.\n\n\n\n\n\ntreewidth(matrix::AbstractMatrix;\n    alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)\n\nA non-mutating version of treewidth!.\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.separator","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.separator","text":"separator(bag::Bag)\n\nGet the separator of a bag.\n\n\n\n\n\nseparator(tree::JunctionTree, i::Integer)\n\nGet the separator at node i.\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.residual","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.residual","text":"residual(bag::Bag)\n\nGet the residual of a bag.\n\n\n\n\n\nresidual(tree::JunctionTree, i::Integer)\n\nGet the residual at node i.\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.relative","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.relative","text":"relative(tree::JunctionTree, i::Integer)\n\nGet the indices in tree[parentindex(tree, i)] corresponding to the elements of separator(tree, i).\n\njulia> using AbstractTrees\n\njulia> using StructuredDecompositions\n\njulia> graph = [\n           0 1 1 0 0 0 0 0\n           1 0 1 0 0 1 0 0\n           1 1 0 1 1 0 0 0\n           0 0 1 0 1 0 0 0\n           0 0 1 1 0 0 1 1\n           0 1 0 0 0 0 1 0\n           0 0 0 0 1 1 0 1\n           0 0 0 0 1 0 1 0\n       ];\n\njulia> label, tree = junctiontree(graph);\n\njulia> bag = tree[parentindex(tree, 1)]\n3-element Bag:\n 6\n 7\n 8\n\njulia> sep = separator(tree, 1)\n2-element view(::Vector{Int64}, 1:2) with eltype Int64:\n 6\n 7\n\njulia> rel = relative(tree, 1)\n2-element view(::Vector{Int64}, 1:2) with eltype Int64:\n 1\n 2\n\njulia> bag[rel] == sep\ntrue\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#Chordal-Graphs","page":"JunctionTrees","title":"Chordal Graphs","text":"","category":"section"},{"location":"api/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"chordalgraph\nischordal\nisperfect","category":"page"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.chordalgraph","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.chordalgraph","text":"chordalgraph([element=true,] tree::JunctionTree)\n\nSee below. The function returns a sparse matrix whose structural nonzeros are filled with element.\n\n\n\n\n\nchordalgraph(Element::Type, tree::JunctionTree)\n\nConstruct the intersection graph of the subtrees of a junction tree. The function returns a sparse matrix with elements of type Element and the same sparsity structure as the lower triangular part of the graph's adjacency matrix.\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.ischordal","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.ischordal","text":"ischordal(matrix::AbstractMatrix)\n\nDetermine whether a simple graph is chordal.\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.isperfect","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.isperfect","text":"isperfect(matrix::AbstractMatrix, order::AbstractVector[, index::AbstractVector])\n\nDetermine whether an fill-reducing permutation is perfect.\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#Elimination-Algorithms","page":"JunctionTrees","title":"Elimination Algorithms","text":"","category":"section"},{"location":"api/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"EliminationAlgorithm\nPermutationOrAlgorithm\nMCS\nRCM\nAMD\nSymAMD\nMMD\nNodeND\nFlowCutter\nSpectral\nBT\npermutation   ","category":"page"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.EliminationAlgorithm","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.EliminationAlgorithm","text":"EliminationAlgorithm\n\nA graph elimination algorithm. The options are\n\ntype name complexity\nMCS maximum cardinality search O(m + n)\nRCM reverse Cuthill-Mckee O(mΔ)\nAMD approximate minimum degree O(mn)\nSymAMD column approximate minimum degree O(mn)\nMMD multiple minimum degree O(mn²)\nNodeND nested dissection \nFlowCutter FlowCutter \nSpectral spectral ordering \nBT Bouchitte-Todinca O(2.6183ⁿ)\n\nfor a graph with m edges, n vertices, and maximum degree Δ.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.PermutationOrAlgorithm","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.PermutationOrAlgorithm","text":"PermutationOrAlgorithm = Union{AbstractVector, EliminationAlgorithm}\n\nEither a permutation or an algorithm.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.MCS","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.MCS","text":"MCS <: EliminationAlgorithm\n\nMCS()\n\nThe maximum cardinality search algorithm.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.RCM","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.RCM","text":"RCM <: EliminationAlgorithm\n\nRCM(; sortbydeg=true)\n\nThe reverse Cuthill-McKee algorithm. Uses SymRCM.jl.\n\nsortbydeg: whether to sort neighbor lists by degree\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.AMD","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.AMD","text":"AMD <: EliminationAlgorithm\n\nAMD(; dense=nothing, aggressive=nothing)\n\nThe approximate minimum degree algorithm. Uses AMD.jl.\n\ndense: dense row parameter\naggressive: aggressive absorbtion\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.SymAMD","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.SymAMD","text":"SymAMD{Index} <: EliminationAlgorithm\n\nSymAMD{Index}(; dense_row=nothing, dense_col=nothing, aggressive=nothing) where Index\n\nSymAMD(; dense_row=nothing, dense_col=nothing, aggressive=nothing)\n\nThe column approximate minimum degree algorithm. Uses AMD.jl.\n\nIndex: either Int or Cint\ndense_row: dense row parameter\ndense_column: dense column parameter\naggressive: aggressive absorbtion\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.MMD","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.MMD","text":"MMD <: EliminationAlgorithm\n\nMMD()\n\nThe multiple minimum degree algorithm. Uses Sparspak.jl.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.NodeND","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.NodeND","text":"NodeND <: EliminationAlgorithm\n\nNodeND()\n\nThe nested dissection algorithm. Uses Metis.jl.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.FlowCutter","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.FlowCutter","text":"FlowCutter <: EliminationAlgorithm\n\nFlowCutter(; time=10, seed=0)\n\nThe FlowCutter algorithm. Uses FlowCutterPACE17_jll.jl. \n\ntime: running time\nseed: random seed\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.Spectral","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.Spectral","text":"Spectral <: EliminationAlgorithm\n\nSpectral(; tol=0.0)\n\nThe spectral ordering algorithm. Uses Laplacians.jl.\n\ntol: tolerance for convergence\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.BT","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.BT","text":"BT <: EliminationAlgorithm\n\nBT()\n\nThe Bouchitte-Todinca algorithm. Uses TreeWidthSolver.jl.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.permutation","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.permutation","text":"permutation(matrix::AbstractMatrix;\n    alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)\n\nConstruct a fill-reducing permutation of the vertices of a simple graph.\n\njulia> using StructuredDecompositions\n\njulia> graph = [\n    0 1 1 0 0 0 0 0\n    1 0 1 0 0 1 0 0\n    1 1 0 1 1 0 0 0\n    0 0 1 0 1 0 0 0\n    0 0 1 1 0 0 1 1\n    0 1 0 0 0 0 1 0\n    0 0 0 0 1 1 0 1\n    0 0 0 0 1 0 1 0\n];\n\njulia> order, index = permutation(graph; alg=MCS());\n\njulia> order\n8-element Vector{Int64}:\n 1\n 6\n 2\n 3\n 4\n 5\n 7\n 8\n\njulia> index == invperm(order)\ntrue\n\n\n\n\n\n","category":"function"},{"location":"api/JunctionTrees/#Supernodes","page":"JunctionTrees","title":"Supernodes","text":"","category":"section"},{"location":"api/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"SupernodeType\nNodal\nMaximal\nFundamental","category":"page"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.SupernodeType","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.SupernodeType","text":"SupernodeType\n\nA type of supernode partition. The options are\n\ntype name\nNodal nodal supernode partition\nMaximal maximal supernode partition\nFundamental fundamental supernode partition\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.Nodal","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.Nodal","text":"Nodal <: SupernodeType\n\nA nodal  supernode partition.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.Maximal","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.Maximal","text":"Maximal <: SupernodeType\n\nA maximal supernode partition.\n\n\n\n\n\n","category":"type"},{"location":"api/JunctionTrees/#StructuredDecompositions.JunctionTrees.Fundamental","page":"JunctionTrees","title":"StructuredDecompositions.JunctionTrees.Fundamental","text":"Fundamental <: SupernodeType\n\nA fundamental supernode partition.\n\n\n\n\n\n","category":"type"},{"location":"api/DecidingSheaves/#DecidingSheaves","page":"DecidingSheaves","title":"DecidingSheaves","text":"","category":"section"},{"location":"api/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"Modules = [StructuredDecompositions.DecidingSheaves]","category":"page"},{"location":"api/DecidingSheaves/#StructuredDecompositions.DecidingSheaves.decide_sheaf_tree_shape","page":"DecidingSheaves","title":"StructuredDecompositions.DecidingSheaves.decide_sheaf_tree_shape","text":"Solve the decision problem encoded by a sheaf.  The algorithm is as follows:    compute on each bag    (optionally, if the decomposition of the solution space   is already known, then it can be passed as an argument),   compute composites on edges,    project back down to bags   answer (providing a witness)     \"no\" if there is an empty bag;      \"yes\" otherwise.\n\n\n\n\n\n","category":"function"},{"location":"api/DecidingSheaves/#StructuredDecompositions.DecidingSheaves.old_adhesion_filter-Tuple{Tuple, StructuredDecompositions.Decompositions.StructuredDecomposition}","page":"DecidingSheaves","title":"StructuredDecompositions.DecidingSheaves.old_adhesion_filter","text":"Filtering algorithm.  Note: we are assuming that we only know how to work with FinSet(Int) !\n\nINPUT: a Finset^{op}-valued structured decomposition d : FG → Span Finset^{op}            (which is expected to be in co-decomposition form;            i.e. as a diagram d : FG → Cospan Finset )           and an indexed span ( (ℓ, r), ( d(ℓ), d(r) ) ) in d            (i.e a pair consisting of span (ℓ, r) in ∫G and its image under d)\n\nOUTPUT: a structured decomposition obtained by replacing the span de in d          by the span obtained by projecting the pullback of de (i.e. taking images)\n\n\n\n\n\n","category":"method"},{"location":"api/DecidingSheaves/#StructuredDecompositions.DecidingSheaves.old_decide_sheaf_tree_shape","page":"DecidingSheaves","title":"StructuredDecompositions.DecidingSheaves.old_decide_sheaf_tree_shape","text":"Solve the decision problem encoded by a sheaf.  The algorithm is as follows:    compute on each bag (optionally, if the decomposition of the solution space                         is already known, then it can be passed as an argument),   compute composites on edges,    project back down to bags   answer (providing a witness)     \"no\" if there is an empty bag;      \"yes\" otherwise.\n\n\n\n\n\n","category":"function"},{"location":"api/Decompositions/#Decompositions","page":"Decompositions","title":"Decompositions","text":"","category":"section"},{"location":"api/Decompositions/","page":"Decompositions","title":"Decompositions","text":"Modules = [StructuredDecompositions.Decompositions]","category":"page"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.StrDecomp","page":"Decompositions","title":"StructuredDecompositions.Decompositions.StrDecomp","text":"Structured decomposition struct     – Consider these as graphs whose vertices are labeled by the objects of some category      and whose edges are labeled by spans in this category\n\n\n\n\n\n","category":"type"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.StrDecomp-Tuple{Any, Any}","page":"Decompositions","title":"StructuredDecompositions.Decompositions.StrDecomp","text":"One can construct a structured decomposition by simply providing the graph representing the shape of the decompostion and the relevant diagram. This constructor will default to setting the Decompsotion Type to Decomposition (i.e. we default to viewing C-valued structured decompositions as diagrams into Span C)\n\n\n\n\n\n","category":"method"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.StrDecomp-Tuple{Catlab.Graphs.BasicGraphs.HasGraph}","page":"Decompositions","title":"StructuredDecompositions.Decompositions.StrDecomp","text":"StrDecomp(graph::HasGraph;\n    alg::PermutationOrAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,\n    snd::SupernodeType=DEFAULT_SUPERNODE_TYPE)\n\nConstruct a structured decomposition of a simple graph. See junctiontree! for the meaning of alg and snd.\n\n\n\n\n\n","category":"method"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.StructuredDecomposition","page":"Decompositions","title":"StructuredDecompositions.Decompositions.StructuredDecomposition","text":"Structured decompositions are diagrams.\n\n\n\n\n\n","category":"type"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.adhesionSpans-Tuple{Any, Any}","page":"Decompositions","title":"StructuredDecompositions.Decompositions.adhesionSpans","text":"Get a vector of indexed adhesion spans; i.e. a vector consisting of pairs (x₁ <- e -> x₂, dx₁ <- de -> dx₂) where x₁ <- e -> x₂ is span in ∫G and dx₁ <- de -> dx₂ is what its image under the decompositon d\n\n\n\n\n\n","category":"method"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.adhesionSpans-Tuple{Any}","page":"Decompositions","title":"StructuredDecompositions.Decompositions.adhesionSpans","text":"Get a vector of the adhesion spans of a decomposition\n\n\n\n\n\n","category":"method"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.adhesions-Tuple{Any, Any}","page":"Decompositions","title":"StructuredDecompositions.Decompositions.adhesions","text":"Get a vector of indexed adhesions; i.e. a vector consisting of pairs (e, de) where e is an edge in ∫G and de is the value mapped to e under the decompositon d\n\n\n\n\n\n","category":"method"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.adhesions-Tuple{Any}","page":"Decompositions","title":"StructuredDecompositions.Decompositions.adhesions","text":"Get a vector of the adhesions of a decomposition\n\n\n\n\n\n","category":"method"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.bags-Tuple{Any, Any}","page":"Decompositions","title":"StructuredDecompositions.Decompositions.bags","text":"Get a vector of indexed bags; i.e. a vector consisting of pairs (x, dx) where x ∈ FG and dx is the value mapped to x under the decompositon d\n\n\n\n\n\n","category":"method"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.bags-Tuple{Any}","page":"Decompositions","title":"StructuredDecompositions.Decompositions.bags","text":"Get a vector of the bags of a decomposition\n\n\n\n\n\n","category":"method"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.∫-Tuple{T} where T<:ACSets.ACSetInterface.ACSet","page":"Decompositions","title":"StructuredDecompositions.Decompositions.∫","text":"Syntactic sugar for costructing the category of elements of a graph.  Note that ∫(G) has type Category whereas elements(G) has type Elements\n\n\n\n\n\n","category":"method"},{"location":"api/Decompositions/#StructuredDecompositions.Decompositions.𝐃","page":"Decompositions","title":"StructuredDecompositions.Decompositions.𝐃","text":"The construction of categories of structured decompostions is functorial;  it consists of a functor 𝐃: Cat{pullback} → Cat taking any category C with pullbacks to the category  𝐃C of C-values structured decompositions. The functoriality of this construction allows us to lift any functor  F : C → E to a functor 𝐃f : 𝐃 C → 𝐃 E which maps C-valued structured decompositions to E-valued structured decompositions.  When we think of the functor F as a computational problem (taking inputs in C to solution spaces in E), then 𝐃f should be thought of as lifting the global comuputation of F to local computation on the constituent parts of C-valued decompositions.  In particular, given a structured decomposition d: FG → C and a sheaf F: C → FinSet^{op} w.r.t to the decompositon topology,  we can make a structured decomposition valued in FinSet^{op} by lifting the sheaf to a functor 𝐃f: 𝐃C → 𝐃(S^{op}) between  categories of structured decompositions. \n\n\n\n\n\n","category":"function"},{"location":"pages/Decompositions/#Decompositions","page":"Decompositions","title":"Decompositions","text":"","category":"section"},{"location":"pages/Decompositions/","page":"Decompositions","title":"Decompositions","text":"Structured decompositions are diagrams. They can be thought of as graphs whose vertices are labeled by the objects of some category and whose edges are labeld by spans in this category. ","category":"page"},{"location":"pages/Decompositions/","page":"Decompositions","title":"Decompositions","text":"abstract type StructuredDecomposition{G, C, D} <: Diagram{id, C, D} end\n\nstruct StrDecomp{G, C, D} <: StructuredDecomposition{G, C, D}  \n  decomp_shape ::G \n  diagram      ::D             \n  decomp_type  ::DecompType\n  domain       ::C  \nend","category":"page"},{"location":"pages/Decompositions/","page":"Decompositions","title":"Decompositions","text":"A structured decomposition can be constructed by providing the graph representing the shap of the decomposition and the relevant diagram. The default constructor will set the decomposition type to Decomposition (viewing C-valued structured decompositions as diagrams into Span C).","category":"page"},{"location":"pages/Decompositions/","page":"Decompositions","title":"Decompositions","text":"StrDecomp(the_decomp_shape, the_diagram)=StrDecomp(the_decomp_shape, the_diagram, Decomposition)","category":"page"},{"location":"pages/Decompositions/","page":"Decompositions","title":"Decompositions","text":"The function StrDecomp will construct a structured decomposition and check whether the decomposition shape makes sense.","category":"page"},{"location":"pages/Decompositions/","page":"Decompositions","title":"Decompositions","text":"StrDecomp(the_decomp_shape, the_diagram, the_decomp_type)","category":"page"},{"location":"pages/Decompositions/","page":"Decompositions","title":"Decompositions","text":"The functions colimit and limit when called on a structured decomposition will take the colimit and limit of the diagram, respectively.","category":"page"},{"location":"pages/Decompositions/","page":"Decompositions","title":"Decompositions","text":"function colimit(d::StructuredDecomposition) colimit(FreeDiagram(d.diagram)) end\nfunction limit(d::StructuredDecomposition) limit(FreeDiagram(d.diagram)) end","category":"page"},{"location":"pages/Decompositions/","page":"Decompositions","title":"Decompositions","text":"The construction of categories of structured decompositions consists of a functor 𝐃:Cat_{pullback} → Cat taking any category C to the category 𝐃C of C-valued structured decompositions. This allows the lifting of any functor F: C → E to a functor 𝐃f : 𝐃C → 𝐃E which maps C-valued structured decompositions to E-valued structured decompostions.","category":"page"},{"location":"pages/Decompositions/","page":"Decompositions","title":"Decompositions","text":"function 𝐃(f, d ::StructuredDecomposition, t::DecompType = d.decomp_type)::StructuredDecomposition","category":"page"},{"location":"pages/JunctionTrees/#JunctionTrees","page":"JunctionTrees","title":"Junction Trees","text":"","category":"section"},{"location":"pages/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"JunctionTrees.jl is a Julia package for constructing tree decompositions of simple graphs. You can use it as follows.","category":"page"},{"location":"pages/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"julia> using StructuredDecompositions\n\njulia> graph = [\n           0 1 1 0 0 0 0 0\n           1 0 1 0 0 1 0 0\n           1 1 0 1 1 0 0 0\n           0 0 1 0 1 0 0 0\n           0 0 1 1 0 0 1 1\n           0 1 0 0 0 0 1 0\n           0 0 0 0 1 1 0 1\n           0 0 0 0 1 0 1 0\n       ];\n\njulia> label, tree = junctiontree(graph);\n\njulia> tree\n6-element JunctionTree:\n[6, 7, 8]\n├─ [1, 6, 7]\n├─ [4, 6, 8]\n│  └─ [3, 4, 6]\n│     └─ [2, 3, 6]\n└─ [5, 7, 8]","category":"page"},{"location":"pages/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"A junction tree is vector of bags, so you can retrieve the bag at node 3 by typing tree[3].","category":"page"},{"location":"pages/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"julia> bag = tree[3]\n3-element Bag:\n 3\n 4\n 6","category":"page"},{"location":"pages/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"Notice that the bag is sorted. Its elements correspond to the vertices label[bag].","category":"page"},{"location":"pages/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"julia> vertices = label[bag]\n3-element Vector{Int64}:\n 7\n 6\n 5","category":"page"},{"location":"pages/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"The width of a junction tree is computed by the function treewidth.","category":"page"},{"location":"pages/JunctionTrees/","page":"JunctionTrees","title":"JunctionTrees","text":"julia> treewidth(tree)\n2","category":"page"},{"location":"pages/DecidingSheaves/#DecidingSheaves","page":"DecidingSheaves","title":"Deciding Sheaves","text":"","category":"section"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"There are two functions that are used here. ","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"The first one called adhesion_filter will take an input Finset^{op}-valued structured decomposition and return a new structured decompostion replacing the span de in d by the span obtained by projecting the pullback of de. ","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"function adhesion_filter(tup::Tuple, d::StructuredDecomposition)","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"The second function is called decidesheaftree_shape and solves a decision problem that is encoded by a sheaf. This function works by computing first on each of the bags, then computes composites on edges and projects back down to bags.","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"function decide_sheaf_tree_shape(f, d::StructuredDecomposition, solution_space_decomp::StructuredDecomposition = 𝐃(f, d, CoDecomposition))","category":"page"},{"location":"pages/DecidingSheaves/#Graph-Coloring-Structure-and-Examples","page":"DecidingSheaves","title":"Graph Coloring Structure and Examples","text":"","category":"section"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"One of the many decision problems that can be solved using this system of functions is figuring out graph colorings. To put it simply, a graph coloring is an assignment of labels(colors) to each vertex of a graph so that no two adjacent vertices have the same label(color).","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"We first define the structure of a coloring.","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"#an H-coloring is a hom onto H\nstruct Coloring\n  n     #the target graph\n  func  #the function mapping opens to lists of homs from G to K_n\nend","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"We then define how we are going to create and test colorings.","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"#construct an n-coloring\nK(n) = complete_graph(Graph, n)\nColoring(n) = Coloring(n, g -> homomorphisms(g, K(n) ))\n#make it callable\n(c::Coloring)(X::Graph) = FinSet(c.func(X)) # given graph homos #f: G₁ → G₂ get morphism col(G₂) → col(G₁) by precomposition: take each λ₂ ∈ col(G₂) to hf ∈ col(G)\nfunction (c::Coloring)(f::ACSetTransformation)  \n  (G₁, G₂)   = (dom(f), codom(f)) \n  (cG₁, cG₂) = (c(G₁), c(G₂))\n  FinFunction( λ₂ -> compose(f,λ₂), cG₂, cG₁ ) #note the contravariance\nend\n\nskeletalColoring(n) = skeleton ∘ Coloring(n)\n\ncolorability_test(n, the_test_case) = is_homomorphic(ob(colimit(the_test_case)), K(n)) == decide_sheaf_tree_shape(skeletalColoring(n), the_test_case)[1]","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"We can consider the example of a ring with seven nodes as our graph.  We first seperate the nodes into bags with adhesions and what our adhesions look like.","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"#bag 1\nH₁ = @acset Graph begin\n    V = 3\n    E = 2\n    src = [1, 2]\n    tgt = [2, 3]\nend\n\n#adhesion 1,2\nH₁₂ = @acset Graph begin\n    V = 2\nend\n\n#bag 2\nH₂ = @acset Graph begin\n    V = 4\n    E = 3\n    src = [1, 2, 3]\n    tgt = [2, 3, 4]\nend\n\nGₛ = @acset Graph begin\n    V = 2\n    E = 1\n    src = [1]\n    tgt = [2]\nend","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"We can then construct our structured decomposition through transformations of these bags and adhesions.","category":"page"},{"location":"pages/DecidingSheaves/","page":"DecidingSheaves","title":"DecidingSheaves","text":"#transformations\nΓₛ⁰ = Dict(1 => H₁, 2 => H₂, 3 => H₁₂)\nΓₛ = FinDomFunctor(\n    Γₛ⁰,\n    Dict(\n      1 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[1], V=[1, 3]),\n      2 => ACSetTransformation(Γₛ⁰[3], Γₛ⁰[2], V=[4, 1]),\n    ),\n    ∫(Gₛ)\n)\n\nmy_decomp1  = StrDecomp(Gₛ, Γₛ)","category":"page"},{"location":"#StructuredDecompositions.jl","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"","category":"section"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"CurrentModule = StructuredDecompositions","category":"page"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"Structural graph theorists and algorithmicists alike know that it's usually a smart idea to decompose graphs into smaller and simpler parts before trying answer difficult computational questions. Tree decompositions are one of the best-known ways of systematically chopping graphs up and, as such they have been key tools for establishing deep results in many areas of discrete mathematics including graph minor theory and algorithmic meta-theorems. ","category":"page"},{"location":"#Generality","page":"StructuredDecompositions.jl","title":"Generality","text":"","category":"section"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"This project consists of an implementation of the category-theoretic notion of structured decompositions. These provide a formalism for decomposing arbitrary mathematical objects (not just graphs) and morphisms between them into smaller constituent parts. Since the definition of a structured decompositions is functorial, one can easily lift computational problems (defined as functors mapping inputs to solution spaces) to functors between from decompositions of the inputs to decompositions of solution spaces.","category":"page"},{"location":"#What-is-in-this-package?","page":"StructuredDecompositions.jl","title":"What is in this package?","text":"","category":"section"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"This package allows one to leverage insights to solve decision problems that are encoded as sheaves efficiently (i.e. in fixed-parameter-tractable time parameterized by the width of the decompositions). Currently, this packages includes many general tools that can be used to decompose arbitrary mathematical objects. One of the many applications of this package, exampled in the Deciding Sheaves module, is graph colorings.","category":"page"},{"location":"","page":"StructuredDecompositions.jl","title":"StructuredDecompositions.jl","text":"Pages = [\n    \"pages/decompostions.md\",\n    \"pages/decidingsheaves.md\",\n    \"pages/junction_trees.md\",\n    \"pages/functor_utils.md\",\n    ]\nDepth = 2","category":"page"}]
}
