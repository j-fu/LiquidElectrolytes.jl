#=
```@meta
Draft=true
```


# ORR cell
([source code](SOURCE_URL))

I-V sweep for Oxygen Reduction

Methods called:
- [`ElectrolyteData`](@@ref)
- [`ivsweep`](@@ref)
- [`dlcapsweep`](@@ref)
- [`PNPSystem`](@@ref)

=#

module Example120_ORRCell
using ExtendableGrids,GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot,Colors, Parameters
using StaticArrays
using LessUnitful
using DocStringExtensions




function main(;voltages=-1:0.005:1,compare=false,molarity=0.1,nref=0,κ=10.0,eneutral=false,scheme=:μex,Plotter=PyPlot,kwargs...)
    
    @local_phconstants R N_A e
    @local_unitfactors nm cm μF mol dm s
    F=N_A*e
    
    
    defaults=(; max_round=3,tol_round=1.0e-10, verbose=false, tol_relative=1.0e-7,tol_mono=1.0e-10)
    kwargs=merge(defaults, kwargs) 
    
    hmin=1.0e-1*nm*2.0^(-nref)
    hmax=1.0*nm*2.0^(-nref)
    L=20.0*nm
    X=geomspace(0,L,hmin,hmax)
    grid=simplexgrid(X)
    

    R0=10.0e-4*mol/(cm^2*s)
    Δg=0.0
    β=0.5
    ϕ_we=0.0
    ihplus=1
    iso4=2
    io2=3

    function halfcellbc(f,u,bnode,data)
        bulkbcondition(f,u,bnode,data)
        @unpack iϕ,eneutral,ϕ_we,Γ_we,RT=data
        if bnode.region==Γ_we
            f.=0.0
            if !data.eneutral
                boundary_dirichlet!(f,u,bnode;species=iϕ,region=data.Γ_we,value=data.ϕ_we)
            end
            μh2o,μ=chemical_potentials!(MVector{4,eltype(u)}(undef),u,data)
            A=(4*μ[ihplus]+μ[io2]-2μh2o+Δg + eneutral*F*(u[iϕ]-data.ϕ_we))/(RT)
            r=rrate(R0,β,A)
            f[ihplus]-=4*r
            f[io2]-=r
        end
    end


    
    
    celldata=ElectrolyteData(;nc=3, z=[1,-2,0], κ=fill(κ,3), Γ_we=1, Γ_bulk=2,eneutral,scheme)

    @unpack iϕ,c_bulk=celldata

    
    c_bulk[io2]=0.001*mol/dm^3
    c_bulk[iso4]=molarity*mol/dm^3
    c_bulk[ihplus]=2.0*molarity*mol/dm^3


    @assert isapprox(celldata.c_bulk'*celldata.z,0, atol=1.0e-12)
    
    cell=PNPSystem(grid;bcondition=halfcellbc,celldata)
    check_allocs!(cell,false)


    ## Compare electroneutral and double layer cases
    if compare
        celldata.eneutral=false
	volts,currs, sols=ivsweep(cell;voltages,ispec=io2,kwargs...)

        celldata.eneutral=true
        nvolts,ncurrs, nsols=ivsweep(cell;voltages,ispec=io2,kwargs...)

        vis=GridVisualizer(;Plotter,resolution=(600,400),clear=true,legend=:lt,xlabel="Δϕ/V",ylabel="I/(A/m^2)")
        scalarplot!(vis,volts,currs,color="red",markershape=:utriangle,markersize=7, markevery=10,label="PNP")
        scalarplot!(vis,nvolts,ncurrs,clear=false,color=:green,markershape=:none,label="NNP")
        return reveal(vis)
    end

    vis=GridVisualizer(resolution=(1200,400),layout=(1,5),Plotter=PyPlot)

    volts,currs, sols=ivsweep(cell;voltages,ispec=io2,kwargs...)
    tsol=VoronoiFVM.TransientSolution(sols,volts)

    for it=1:length(tsol.t)
        tsol.u[it][io2,:]/=mol/dm^3
        tsol.u[it][ihplus,:]/=mol/dm^3
        tsol.u[it][iso4,:]/=mol/dm^3
    end

    xmax=20*nm
    xlimits=[0,xmax]
    aspect=3.5*xmax/(tsol.t[end]-tsol.t[begin])

    scalarplot!(vis[1,1],currs,volts,markershape=:none,title="IV",xlabel="I",ylabel="V")
    scalarplot!(vis[1,2],cell,tsol;species=io2,aspect,xlimits,title="O2",colormap=:summer)
    scalarplot!(vis[1,3],cell,tsol;species=ihplus,aspect,xlimits,title="H+",colormap=:summer)
    scalarplot!(vis[1,4],cell,tsol;species=iϕ,aspect,xlimits,title="ϕ",colormap=:bwr)
    scalarplot!(vis[1,5],tsol[io2,1,:],volts,title="c_o2(0)",xlabel="O2",ylabel="V")

    reveal(vis)
end

end
