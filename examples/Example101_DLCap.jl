#=

# Double Layer Capacitance

([source code](@__SOURCE_URL__))

Calculation of double layer capacitance of a symmetric 1:1 electrolyte.

![](Example101_DLCap.svg)

Methods called:
- [`ElectrolyteData`](@@ref)
- [`dlcapsweep`](@@ref)
- [`PNPSystem`](@@ref)
=#
module Example101_DLCap
using DynamicQuantities
using VoronoiFVM,ExtendableGrids,GridVisualize
using LiquidElectrolytes
using Colors


function main(;voltages=collect(-2:0.01:2)*u"V",           ## Voltages/V
              molarities=[0.001,0.01,0.1,1]*u"mol/dm^3", ## Molarities/M
	      nref=0,	                     ## Refinement level
	      scheme=:μex,	             ## Flux calculation scheme
	      κ=10.0,                        ## Solvation number
              Plotter=nothing,               ## Plotter
	      kwargs...                      ## Solver kwargs
              )
    
    ## Obtain unit factors from LessUnitful.jl

    nm=u"nm"|>ustrip
    cm=u"cm"|>ustrip
    m=u"m"|>ustrip
    μF=u"μF"|>ustrip
    V=u"V"|>ustrip
    M=u"mol/dm^3"|>ustrip
    
    ## Create a standard 1D grid with grid spacing following geometric progression
    hmin=1.0e-1*nm*2.0^(-nref)
    hmax=1.0*nm*2.0^(-nref)
    L=20.0*nm
    X=geomspace(0,L,hmin,hmax)
    grid=simplexgrid(X)

    ## Define boundary conditions 
    function bcondition(f,u,bnode,data)
	(;iϕ,Γ_we,Γ_bulk,ϕ_we) = data
        
	## Dirichlet ϕ=ϕ_we at Γ_we
	boundary_dirichlet!(f,u,bnode,species=iϕ,region=Γ_we,value=ϕ_we)

        ## Bulk condition at Γ_bulk
        bulkbcondition(f,u,bnode,data,region=Γ_bulk)
    end

    ## Create electrolyte data
    celldata=ElectrolyteData(nc=2,
                             Γ_we=1,
			     Γ_bulk=2;
			     scheme,
                             κ=fill(κ,2),
                             c_bulk=fill(0.01,2)*u"mol/dm^3")
    
    ## Create Poisson-Nernst-Planck system
    cell=PNPSystem(grid;bcondition,celldata)

    @show electrolytedata(cell)
    
    ## Visualization
    vis=GridVisualizer(;resolution=(500,300),
                       legend=:rt,
                       clear=true,
                       xlabel="φ/V",
                       ylabel="C_dl/(μF/cm^2)",
                       Plotter)

    for imol=1:length(molarities)
        
	color=RGB(1-imol/length(molarities),0,imol/length(molarities))

        
	result=dlcapsweep(cell;
                          δ=1.0e-6,
                          voltages=collect(voltages),
                          molarity=molarities[imol],
                          kwargs...)

	cdl0=dlcap0(celldata)
	scalarplot!(vis,
                    result.voltages .|> ustrip,
                    result.dlcaps .|> uconvert(us"μF/cm^2") .|>ustrip;
                    color,
                    clear=false,
                    label="$(molarities[imol])M",
                    markershape=:none)
        
	scalarplot!(vis,[0],[cdl0].|>uconvert(us"μF/cm^2")|>ustrip;
                    clear=false,
                    markershape=:circle,
                    label="")
    end
    reveal(vis)
end

function generateplots(dir; Plotter = nothing, kwargs...)    #hide
    if ismakie(Plotter)                                      #hide
        Plotter.activate!(; type = "svg", visible = false)   #hide
        p=main(;Plotter)                                     #hide
        Plotter.save(joinpath(dir, "Example101_DLCap.svg"), p)  #hide
    end                                                      #hide
    nothing                                                  #hide
end                                                          #hide
end 
