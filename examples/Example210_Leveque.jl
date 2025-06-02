module Example210_Leveque

using VoronoiFVM
using Printf
using ExtendableGrids
using DelimitedFiles
using GridVisualize
using LessUnitful
using Test

function create_peak_velocity(coordsystem = Cartesian2D; LPe = 8.0, Hcell = 1.0, Diff = 1.0)
    if coordsystem <: Cartesian2D    
        return Pe -> Pe*Hcell*Diff/(4.0*LPe^2)
    elseif coordsystem <: Cylindrical2D
        return Pe -> (Pe * Diff) / (2 * Hcell)
    end
end

function create_velocity_field(coordsystem = Cartesian2D; Hcell = 1.0)
    if coordsystem <: Cartesian2D
        return (x,y) -> (4*(y/Hcell)*(1.0-(y/Hcell)),0)
    elseif coordsystem <: Cylindrical2D
        return (x,y) -> (0.0, 2 * ((Hcell)^2 - x^2)/((Hcell)^2))
    end
end

function create_grid(coordsystem = Cartesian2D, Hcell = 1.0, nref = 0, Lbreak = 1.0, Lcat = 8.0; binert = 1, bin = 3, bout = 2, brea = 4, brot = 5)
    @local_unitfactors nm μm mm m s
    Lcell=2*Lbreak+Lcat  #overall length of cell 
    hy_electrode=0.01*Hcell/2.0^nref
    hy_top=0.2*Hcell/2.0^nref
    hx_electrode=0.01*Lbreak/2.0^nref
    hx_max=0.3*Lbreak/2.0^nref

    xcoord0=geomspace(0.0,Lbreak,hx_max,hx_electrode)
    xcoord1=geomspace(Lbreak,Lcat+Lbreak,hx_electrode,hx_max)
    xcoord2=geomspace(Lcat+Lbreak,Lcat+2*Lbreak,hx_max,hx_max)
    xcoord=glue(xcoord0,xcoord1)
    xcoord=glue(xcoord,xcoord2)
    ycoord=geomspace(0.0,Hcell,hy_electrode,hy_top)
    grid=simplexgrid(xcoord,ycoord)
    cellmask!(grid,[0,0],[Lcell,Hcell],1, tol=1*nm)

    bfacemask!(grid,[0,0],[Lcell,0],binert, tol=1*nm)
    bfacemask!(grid,[0,Hcell],[Lcell,Hcell],binert, tol=1*nm)

    bfacemask!(grid,[0,0],[0,Hcell],bin, tol=1*nm)
    bfacemask!(grid,[Lcell,0],[Lcell,Hcell],bout, tol=1*nm)
    bfacemask!(grid,[Lbreak,0],[Lcell-Lbreak,0],brea, tol=1*nm)

    if coordsystem <: Cylindrical2D
        bfacemask!(grid,[0,Hcell],[Lcell,Hcell],brot, tol=1*nm)
        grid[Coordinates][1,:], grid[Coordinates][2,:] = -grid[Coordinates][2,:], grid[Coordinates][1,:]
        grid[Coordinates][1,:] .+= Hcell
        grid = circular_symmetric!(grid)
    end

    return grid
end

function create_ShFac(coordsystem = Cartesian2D; Diff = 1.0, cin = 1.0, LPe = 8.0)
    if coordsystem <: Cartesian2D
        return (Diff*(cin-0))
    elseif coordsystem <: Cylindrical2D
        # see Newman, 1969: https://doi.org/10.1115/1.3580091
        return π * LPe * Diff * (cin - 0)
    end
end

function create_asymptotics(coordsystem = Cartesian2D; LPe = 8.0, Hcell = 1.0)
    if coordsystem <: Cartesian2D
        return x -> (0.8075491*(x^(1.0/3.0)))
    elseif coordsystem <: Cylindrical2D
        return x -> (1.6151 * ((x / (LPe / (2 * Hcell)))^(1.0 / 3.0)))
    end
end

function main(coordsystem=Cartesian2D;
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

    cylindrical = coordsystem <: Cylindrical2D

    brot = 5
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

    LPe=Lcat # length value for Peclet number calculation

    grid = create_grid(coordsystem, Hcell, nref, Lbreak, Lcat; binert, bin, bout, brea, brot)

    if verbose
        println(grid)
    end
    
    if plotgrid
        gridplot(grid,Plotter=Plotter)
        return
    end

    Velocity = create_peak_velocity(coordsystem; LPe, Hcell, Diff)
    vel = create_velocity_field(coordsystem; Hcell)

    function diffusion(f,u,edge, data)
        f[1]=Diff*(u[1]-u[2])
    end
    peak_velocity::Float64=0.0
    evelo=edgevelocities(grid,vel)
    bfvelo=bfacevelocities(grid,vel)

    function convdiff(f,u,edge, data)
        vd=peak_velocity*evelo[edge.index]/Diff
        bp=fbernoulli(vd)
        bm=fbernoulli(-vd)
        bp=bp*Diff
        bm=bm*Diff
        f[1]=bp*u[1] - bm*u[2]
    end
    
    function outflow(f,u,node, data)
        if node.region==bout
            f[1]=peak_velocity*bfvelo[node.ibnode,node.ibface]*u[1]
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

    ShFac = create_ShFac(coordsystem; Diff, cin, LPe)

    while Pe<PeMax
        peak_velocity=Velocity(Pe)
        solution=solve(sys;inival,control=ctl)
        I=VoronoiFVM.integrate(sys,tfc_rea,solution)
        Iin=VoronoiFVM.integrate(sys,tfc_in,solution)
        Iout=VoronoiFVM.integrate(sys,tfc_out,solution)
        Sh=abs(I[ispec]/ShFac)
        if verbose
            @printf("Pe=%.3e Sh=%.3e v=%.3e I: %.3e %.3e %.3e sum=%.3e\n",Pe,Sh, peak_velocity,I[ispec],Iin[ispec], Iout[ispec],I[ispec]+Iin[ispec]+Iout[ispec] )
        end
        push!(Pes,Pe)
        push!(Shs,Sh)
        if plotsolution
            scalarplot!(vis[1,1],grid,solution[1,:],title=@sprintf("velo=%.2g\n",peak_velocity),aspect=3,show=true)
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

        leveque = create_asymptotics(coordsystem; LPe, Hcell)
        
        vis=GridVisualizer(;Plotter, xscale=:log, yscale=:log, xlabel="Pe", ylabel="Sh", legend=:lt)
        scalarplot!(vis, Pes,leveque.(Pes), color=:darkgreen,label="Leveque asymptotics")
        scalarplot!(vis, Pes,Shs, color=:red, markevery=1, markersize=8,label="RRDE-Julia", clear=false, markershape=:circle)
        if !cylindrical
            scalarplot!(vis, refdata[1,:],refdata[2,:], color=:darkblue,label="pdelib", clear=false, markershape=:none)
        end
        reveal(vis)
        save(pesh,vis)
    end
    return sum(Shs)
end

function runtests()
    @test main(Cartesian2D; microscale=false, verbose=false)≈2.1918830162624043e8
    @test main(Cartesian2D; microscale=true, verbose=false)≈716700.6425206012
    @test main(Cylindrical2D; microscale=false, verbose=false)≈7.001740151039522e9
    @test main(Cylindrical2D; microscale=true, verbose=false)≈3.500919519297249e8
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
