push!(LOAD_PATH,joinpath(@__DIR__,".."))
push!(LOAD_PATH,joinpath(@__DIR__,"..","examples"))
using Documenter, LiquidElectrolytes, LessUnitful,Literate, Revise

function mkdocs()
    example_jl_dir = joinpath(@__DIR__,"..","examples")
    example_md_dir  = joinpath(@__DIR__,"src","examples")

    rm(example_md_dir,force=true,recursive=true)
    
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
    
    DocMeta.setdocmeta!(LiquidElectrolytes, :DocTestSetup, :(using LiquidElectrolytes, Unitful, LessUnitful); recursive=true)

    makedocs(sitename="LiquidElectrolytes.jl",
             modules = [LiquidElectrolytes],
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
             ])
    if !isinteractive()
        deploydocs(repo = "github.com/j-fu/LiquidElectrolytes.jl.git", devbranch = "main")
    end

end

mkdocs()

