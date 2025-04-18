module Example210_Leveque

using VoronoiFVM
using Printf
using ExtendableGrids
using DelimitedFiles
using GridVisualize
using LessUnitful
using Test

function main(;
              nref=0,
              Plotter=nothing,
              plotgrid=false,
              plotsolution=false,
              microscale=false,
              verbose=true,
              xyscale=1.0,
              Pe=1.0,
              PeMax=1.0e15,
              dir="./"
              )
    @local_unitfactors nm μm mm m s
    
    binert=1
    bin=3
    bout=2
    brea=4
    cin=1
    Diff::Float64=0.0
    
    if microscale
        Diff=9.5e-9*m^2/s # diffusion coefficient 
        Lcat=8*mm     # length of catalyst
        Hcell=50*μm        # cell height
        Lbreak=1*mm    # length of boundary before and after catalyst
    else
        Lcat=8*xyscale
        Hcell=1*xyscale
        Diff=1*sqrt(xyscale)
        Lbreak=1*xyscale
    end


    Lcell=2*Lbreak+Lcat  #overall length of cell 
    LPe=Lcat         # length value for Peclet number calculation

    hy_electrode=0.01*Hcell/2.0^nref
    hy_top=0.2*Hcell/2.0^nref
    hx_electrode=0.01*Lbreak/2.0^nref
    hx_max=0.3*Lbreak/2.0^nref
    
    function Velocity(Pe)
        return Pe*Hcell*Diff/(4.0*LPe^2)
    end

    function fhp(x,y)
        yh=y/Hcell
        return 4*yh*(1.0-yh),0
    end

    xcoord0=geomspace(0.0,Lbreak,hx_max,hx_electrode)
    xcoord1=geomspace(Lbreak,Lcat+Lbreak,hx_electrode,hx_max)
    xcoord2=geomspace(Lcat+Lbreak,Lcat+2*Lbreak,hx_max,hx_max)
    xcoord=glue(xcoord0,xcoord1)
    xcoord=glue(xcoord,xcoord2)
    ycoord=geomspace(0.0,Hcell,hy_electrode,hy_top)

    
    grid=simplexgrid(xcoord,ycoord)
    if verbose
        println(grid)
    end
    cellmask!(grid,[0,0],[Lcell,Hcell],1, tol=1*nm)
    bfacemask!(grid,[0,0],[Lcell,Hcell],binert, tol=1*nm)
    bfacemask!(grid,[0,0],[0,Hcell],bin, tol=1*nm)
    bfacemask!(grid,[Lcell,0],[Lcell,Hcell],bout, tol=1*nm)
    bfacemask!(grid,[Lbreak,0],[Lcell-Lbreak,0],brea, tol=1*nm)
    
    if plotgrid
        gridplot(grid,Plotter=Plotter)
        return
    end

    function diffusion(f,u,edge, data)
        f[1]=Diff*(u[1]-u[2])
    end
    Velo::Float64=0.0
    evelo=edgevelocities(grid,fhp)
    bfvelo=bfacevelocities(grid,fhp)

    function convdiff(f,u,edge, data)
        vd=Velo*evelo[edge.index]/Diff
        bp=fbernoulli(vd)
        bm=fbernoulli(-vd)
        bp=bp*Diff
        bm=bm*Diff
        f[1]=bp*u[1] - bm*u[2]
    end
    
    function outflow(f,u,node, data)
        if node.region==bout
            f[1]=Velo*bfvelo[node.ibnode,node.ibface]*u[1]
        end
    end
    
    ispec=1
    physics=VoronoiFVM.Physics(num_species=1,flux=convdiff,breaction=outflow)
    sys=VoronoiFVM.DenseSystem(grid,physics)
    enable_species!(sys,ispec,[1])


    boundary_dirichlet!(sys,ispec,bin,cin)
    boundary_dirichlet!(sys,ispec,brea,0.0)
    
    factory=VoronoiFVM.TestFunctionFactory(sys)
    tfc_rea=testfunction(factory,[bin,bout],[brea])
    tfc_in=testfunction(factory,[bout,brea],[bin])
    tfc_out=testfunction(factory,[bin,brea],[bout])
    
    
    inival=unknowns(sys)
    solution=unknowns(sys)
    inival.=0
    solution.=0
    
    if plotsolution
        vis=GridVisualizer(Plotter=Plotter)
    end

    Pes=zeros(0)
    Shs=zeros(0)
 
    ctl=VoronoiFVM.NewtonControl()

    while Pe<PeMax
        Velo=Velocity(Pe)
        solution=solve(sys;inival,control=ctl)
        ShFac=(Diff*(cin-0))
        I=VoronoiFVM.integrate(sys,tfc_rea,solution)
        Iin=VoronoiFVM.integrate(sys,tfc_in,solution)
        Iout=VoronoiFVM.integrate(sys,tfc_out,solution)
        Sh=abs(I[ispec]/ShFac)
        if verbose
            @printf("Pe=%.3e Sh=%.3e v=%.3e I: %.3e %.3e %.3e sum=%.3e\n",Pe,Sh, Velo,I[ispec],Iin[ispec], Iout[ispec],I[ispec]+Iin[ispec]+Iout[ispec] )
        end
        push!(Pes,Pe)
        push!(Shs,Sh)
        if plotsolution
            scalarplot!(vis[1,1],grid,solution[1,:],title=@sprintf("velo=%.2g\n",Velo),aspect=3,show=true)
        end
        Pe*=2
        inival.=solution
    end

    if !isnothing(Plotter)
        pdelib_data=joinpath(@__DIR__,"..","assets","pdelib-macro.dat")
        pesh=joinpath(dir,"leveque-macro.png")
        if microscale
            pdelib_data=joinpath(@__DIR__,"..","assets","pdelib-micro.dat")
            pesh=joinpath(dir,"leveque-micro.png")
        end
        
        refdata=transpose(readdlm(pdelib_data, comments=true, comment_char='#'))
        function leveque(x)
            return 0.8075491*(x^(1.0/3.0))
        end
        
        vis=GridVisualizer(;Plotter, xscale=:log, yscale=:log, xlabel="Pe", ylabel="Sh", legend=:lt)
        scalarplot!(vis, Pes,leveque.(Pes), color=:darkgreen,label="Leveque asymptotics")
        scalarplot!(vis, Pes,Shs, color=:red, markevery=1, markersize=8,label="RRDE-Julia", clear=false, markershape=:circle)
        scalarplot!(vis, refdata[1,:],refdata[2,:], color=:darkblue,label="pdelib", clear=false, markershape=:none)
        reveal(vis)
        save(pesh,vis)
    end
    return sum(Shs)
end

function runtests()
    @test main(microscale=false, verbose=false)≈2.1918830162624043e8
    @test main(microscale=true, verbose=false)≈716700.6425206012
end

function generateplots(dir; Plotter = nothing, kwargs...)    #hide
    if ismakie(Plotter)                                      #hide
        Plotter.activate!(; type = "png", visible = false)   #hide
        main(;Plotter, microscale=false,dir)          #hide
        main(;Plotter, microscale=true,dir)           #hide
    end                                                      #hide
    nothing                                                  #hide
end                                                          #hide

end
