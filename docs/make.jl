
push!(LOAD_PATH,joinpath(@__DIR__,".."))
push!(LOAD_PATH,joinpath(@__DIR__,"..","examples"))
using Documenter, LiquidElectrolytes, LessUnitful,Literate, Revise, PlutoStaticHTML


const NOTEBOOK_DIR = joinpath(@__DIR__, "..", "notebooks")
const NOTEBOOKS=["DLCap", "ORR"]
const NOTEBOOKS_JL=NOTEBOOKS.*".jl"
const NOTEBOOKS_MD=NOTEBOOKS.*".md"


function build_all_notebooks()
    thisdir=pwd()
    Pkg.activate(NOTEBOOK_DIR)
    Pkg.instantiate()
    Pkg.activate(thisdir)
    println("Building notebooks in $NOTEBOOK_DIR")
    ENV["PLUTO_PROJECT"]=NOTEBOOK_DIR
    oopts = OutputOptions(; append_build_context=true)
    output_format = documenter_output
    bopts = BuildOptions(NOTEBOOK_DIR; output_format)
    build_notebooks(bopts,NOTEBOOKS_JL, oopts)
    return nothing
end

function mkdocs()
    example_jl_dir = joinpath(@__DIR__,"..","examples")
    example_md_dir  = joinpath(@__DIR__,"src","examples")
    notebook_md_dir  = joinpath(@__DIR__,"src","notebooks")

    
    
    rm(example_md_dir,force=true,recursive=true)
    rm(notebook_md_dir,force=true,recursive=true)
    mkdir(notebook_md_dir)
    
    function replace_source_url(input,source_url)
        lines_in = collect(eachline(IOBuffer(input)))
        lines_out=IOBuffer()
        for line in lines_in
            println(lines_out,replace(line,"SOURCE_URL" => source_url))
        end
        return String(take!(lines_out))
    end

    function replace_atat(input)
        lines_in = collect(eachline(IOBuffer(input)))
        lines_out=IOBuffer()
        for line in lines_in
            println(lines_out,replace(line,"@@" => "@"))
        end
        return String(take!(lines_out))
    end
    wd=joinpath(pwd(),"..")
    for example_source in readdir(example_jl_dir)
        base,ext=splitext(example_source)
        if ext==".jl"
            source_url="https://github.com/j-fu/LiquidElectrolytes.jl/blob/main/examples/"*example_source
            preprocess(buffer)=replace_source_url(buffer,source_url)
            Literate.markdown(joinpath(@__DIR__,"..","examples",example_source),
                              example_md_dir;
                              documenter=true,
                              info=false,
                              preprocess,
                              postprocess=replace_atat)
        end
    end
    
    
    #generated_examples=vcat(["runexamples.md"],joinpath.("examples",readdir(example_md_dir)))
    generated_examples=joinpath.("examples",readdir(example_md_dir))
    
    
    build_all_notebooks()
    for nb in NOTEBOOKS_MD
        mv(joinpath(NOTEBOOK_DIR,nb),joinpath(notebook_md_dir,nb))
    end
    notebooks=joinpath.("notebooks",NOTEBOOKS_MD)

    notebooks=[ nb*".jl"=> joinpath("notebooks",nb*".md") for nb in NOTEBOOKS ]
    
    DocMeta.setdocmeta!(LiquidElectrolytes, :DocTestSetup, :(using LiquidElectrolytes, Unitful, LessUnitful); recursive=true)

    makedocs(sitename="LiquidElectrolytes.jl",
             modules = [LiquidElectrolytes],
             format = Documenter.HTML(mathengine=MathJax3()),
             clean = false, 
             doctest = true,
             draft = false,
             authors = "J. Fuhrmann",
             repo="https://github.com/j-fu/LiquidElectrolytes.jl/",
             pages=[
                 "Home"=>"index.md",
                 "Electrolyte models"=>"api.md",
                 "Standard calculations"=>"std.md",
                 "Changes" => "changes.md",
                 "Examples" => generated_examples,
                 "Internal API"=>"internal.md",
                 "Notebooks" => notebooks,
             ])
    if !isinteractive()
        deploydocs(repo = "github.com/j-fu/LiquidElectrolytes.jl.git", devbranch = "main")
    end

end

mkdocs()

