push!(LOAD_PATH, joinpath(@__DIR__, ".."))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "examples"))

using Documenter, ExampleJuggler, CairoMakie, LiquidElectrolytes, PlutoStaticHTML, VoronoiFVM

ExampleJuggler.verbose!(true)
thisdir = pwd()

function make(; with_notebooks = true, with_examples = true)
    pages = Any[
        "index.md",
        "notations.md",
        "api.md",
        "std.md",
        "internal.md",
    ]

    cleanexamples()
    exampledir = joinpath(@__DIR__, "..", "examples")
    notebookdir = joinpath(@__DIR__, "..", "notebooks")


    size_threshold_ignore = []
    if with_notebooks
        notebooks = [
            "EquilibriumCheck.jl",
            "ORR.jl",
            "ElectroOsmosis.jl",
        ] #, "BufferReactions.jl", "SurfaceKinetics_draft.jl"]
        notebook_examples = @docplutonotebooks(notebookdir, notebooks, iframe = false, append_build_context = false)
        size_threshold_ignore = last.(notebook_examples)
        push!(pages, "Notebooks" => notebook_examples)
    end

    if with_examples
        modules = [
            "Example101_DLCap.jl",
            "Example110_Fe23Cell.jl",
            "Example120_ORRCell.jl",
        ]
        module_examples = @docmodules(exampledir, modules, use_module_titles = true, Plotter = CairoMakie)
        push!(pages, "Examples" => module_examples)
    end


    DocMeta.setdocmeta!(LiquidElectrolytes, :DocTestSetup, :(using LiquidElectrolytes, Unitful, LessUnitful); recursive = true)

    cd(thisdir)
    makedocs(;
        sitename = "LiquidElectrolytes.jl",
        modules = [LiquidElectrolytes],
        format = Documenter.HTML(; mathengine = MathJax3(), size_threshold_ignore),
        clean = false,
        doctest = true,
        warnonly = true,
        draft = false,
        authors = "J. Fuhrmann",
        repo = "https://github.com/j-fu/LiquidElectrolytes.jl/",
        pages
    )

    cleanexamples()
    if !isinteractive()
        deploydocs(repo = "github.com/j-fu/LiquidElectrolytes.jl.git", devbranch = "main")
    end
    return
end

if isinteractive() || (haskey(ENV, "DOCSONLY") && ENV["DOCSONLY"] == "true")
    make(; with_notebooks = false, with_examples = false)
else
    make()
end
