module Example211_Levich
using VoronoiFVM
using ExtendableGrids
using Printf
using LinearAlgebra
using GridVisualize
using LiquidElectrolytes.RotatingElectrodes
using LessUnitful
using Test

"""
    main(;
              nref=0,
              plot_grid=false,
              plot_conc=false,
              disktype="E6R1",
              Diff=1.0e-6,
              minfreq=1.0,
              maxfreq=256, 
              Plotter=nothing,
              verbose=true
              )


"""
function main(;
              nref=0,
              plot_grid=false,
              plot_conc=false,
              plot_results=false,
              disktype="E6R1",
              Diff=1.0e-6,
              minfreq=1.0,
              maxfreq=256, 
              Plotter=nothing,
              verbose=true,
              adjust_grid=true
              )
    @local_unitfactors nm μm mm cm m s mA
    FaradayConstant=ph"AvogadroConstant"*ph"ElementaryCharge"
    
    geom=DiscSpec(disktype)
    if adjust_grid
        hz=(0.05,0.05)
    else
        hz=(0.005,0.1)
    end
    grid=rrde_grid(geom;nref, hz)
    grid_points=num_nodes(grid)

    if plot_grid
        gridplot(grid,Plotter=Plotter,edges=false, linewidth=0.1)
        return
    end
    
    
    
    nu=0.02942*cm^2/s
    Diff=Diff*cm^2/s
    cin=1
    rrde_data=RRDEData([Diff])
    
    nfreq1=Int32(ceil(sqrt(maxfreq)))
    nfreq0=Int32(ceil(sqrt(minfreq)))
    
    ispec=1
    karman_data=KarmanData(1.0,1.0)
    physics=VoronoiFVM.Physics(num_species=1,flux=rrde_flux,data=rrde_data,breaction=rrde_outflow)
    rdcell=VoronoiFVM.DenseSystem(grid,physics)
    enable_species!(rdcell,ispec,[1])
    inival=unknowns(rdcell)
    inival.=0
    sol=unknowns(rdcell)
    
    factory=VoronoiFVM.TestFunctionFactory(rdcell)
    tfc_disk=testfunction(factory,[b_in,b_out,b_ring],[b_disk])
    tfc_ring=testfunction(factory,[b_in,b_out,b_disk],[b_ring])
    tfc_in=testfunction(factory,[b_out,b_disk,b_ring],[b_in])
    tfc_out=testfunction(factory,[b_in,b_disk,b_ring],[b_out])
    control=VoronoiFVM.NewtonControl()

    
    freqs=[]
    res_idisk_sim=[]
    res_idisk_levich=[]
    res_idisk_newman=[]
    res_iring_sim=[]
    res_iring_levich=[]
    res_coleff_sim=[]
    res_coleff_albery=[]

    if plot_conc
        vis=GridVisualizer(layout=(1,3),Plotter=Plotter,resolution=(1200,400),fignumber=4)
    end
    
    for ifreq=nfreq0:nfreq1
        freq=ifreq*ifreq
        push!(freqs,freq)
        omega=freq*2.0*π
        @printf("---- rotation frequency: %.2f omega: %.2e  -------\n",freq,omega)

        if adjust_grid
            #!!     adjustment of the grid relative to rotation frequency
            delta=1.61*(nu)^(1/6)*(Diff)^(1/3)/sqrt(omega)
            rescale_z!(grid,delta)
            update_grid!(rdcell)
        end
        xmin=minimum(grid[Coordinates][1,:])
        xmax=maximum(grid[Coordinates][1,:])
        zmin=minimum(grid[Coordinates][2,:])
        zmax=maximum(grid[Coordinates][2,:])
        aspect=0.5*(xmax-xmin)/(zmax-zmin)
        
        KarmanData!(karman_data,omega,nu)
        fkarman(r,z)=fKarman(karman_data,r,z)
        RotatingElectrodes.update!(rrde_data,grid,fkarman)

        ################# DISK LIMITING CURRENT
        # fixed concentration at inlet
        boundary_dirichlet!(rdcell,ispec,b_in,cin)

        #  homogeneous Dirichlet at disk
        boundary_dirichlet!(rdcell,ispec,b_disk,0)

        #  homogeneous Neumann at ring
        boundary_neumann!(rdcell,ispec,b_ring,0)

        state=VoronoiFVM.SystemState(rdcell)
        sol=solve!(state;inival)

        I_disk=VoronoiFVM.integrate(rdcell,tfc_disk,sol)
        I_ring=VoronoiFVM.integrate(rdcell,tfc_ring,sol)
        I_in=VoronoiFVM.integrate(rdcell,tfc_in,sol)
        I_out=VoronoiFVM.integrate(rdcell,tfc_out,sol)
        if verbose
        @printf("[d] ref=%2d f=%3d in=%+8.3e out=%+8.3e disk=%+8.3e ring=%+8.3e bal=%+8.3e max-cin=%+8.3e \n",
                nref,
                freq,
                I_in[1],
                I_out[1],
                I_disk[1],
                I_ring[1],
                I_in[1]+I_out[1]+I_disk[1]+I_ring[1],
                maximum(sol)-cin
                )
        end
        if plot_conc
            vis[1,1][:grid]=nothing
            scalarplot!(vis[1,1],grid,sol[1,:],aspect=aspect,title="idisk")
        end
        push!(res_idisk_sim,abs(FaradayConstant*I_disk[1]))
        push!(res_idisk_levich,idisk_levich(geom,omega,nu,Diff))
        push!(res_idisk_newman,idisk_newman(geom,omega,nu,Diff))

        ################# RING LIMITING CURRENT
        # fixed concentration at inlet
        boundary_dirichlet!(rdcell,ispec,b_in,cin)

        #  homogeneous Neumann at disk
        boundary_neumann!(rdcell,ispec,b_disk,0)

        #  homogeneous Dirichlet at ring
        boundary_dirichlet!(rdcell,ispec,b_ring,0)


        sol=solve!(state;inival)

        I_disk=VoronoiFVM.integrate(rdcell,tfc_disk,sol)
        I_ring=VoronoiFVM.integrate(rdcell,tfc_ring,sol)
        I_in=VoronoiFVM.integrate(rdcell,tfc_in,sol)
        I_out=VoronoiFVM.integrate(rdcell,tfc_out,sol)

        if verbose
        @printf("[r] ref=%2d f=%3d in=%+8.3e out=%+8.3e disk=%+8.3e ring=%+8.3e bal=%+8.3e max-cin=%+8.3e \n",
                nref,
                freq,
                I_in[1],
                I_out[1],
                I_disk[1],
                I_ring[1],
                I_in[1]+I_out[1]+I_disk[1]+I_ring[1],
                maximum(sol)-cin
                )
        end
        if plot_conc
            vis[1,2][:grid]=nothing
            scalarplot!(vis[1,2],grid,sol[1,:],aspect=aspect,title="iring")
        end

        push!(res_iring_sim,abs(FaradayConstant*I_ring[1]))
        push!(res_iring_levich,iring_levich(geom,omega,nu,Diff))

        ################# COLLECTION EFFICIENCY
        # homogeneous Dirichlet at inlet
        # The rationale here is that stuff diffusing away that far
        # would essentially leave the cell anyway.
        # Neumann would force it to stay in the cell, slightly increasing
        # the collection efficieny.
        
        boundary_dirichlet!(rdcell,ispec,b_in,0.0)

        #  inhomogeneous Neumann at disk
        boundary_neumann!(rdcell,ispec,b_disk,1.0e-5)

        #  homogeneous Dirichlet at ring
        boundary_dirichlet!(rdcell,ispec,b_ring,0)
        
        sol=solve!(state;inival)


        I_disk=VoronoiFVM.integrate(rdcell,tfc_disk,sol)
        I_ring=VoronoiFVM.integrate(rdcell,tfc_ring,sol)
        I_in=VoronoiFVM.integrate(rdcell,tfc_in,sol)
        I_out=VoronoiFVM.integrate(rdcell,tfc_out,sol)

        if verbose
        @printf("[c] ref=%2d f=%3d in=%+8.3e out=%+8.3e disk=%+8.3e ring=%+8.3e bal=%+8.3e max-cin=%+8.3e \n",
                nref,
                freq,
                I_in[1],
                I_out[1],
                I_disk[1],
                I_ring[1],
                I_in[1]+I_out[1]+I_disk[1]+I_ring[1],
                maximum(sol)-cin
                )
        end
        if plot_conc
            vis[1,3][:grid]=nothing
            scalarplot!(vis[1,3],grid,sol[1,:],aspect=aspect,title="coleff")
            reveal(vis)
        end
        push!(res_coleff_sim,abs(I_ring[1]/I_disk[1]))
        push!(res_coleff_albery,coleff_albery(geom))
        println("")
        
    end

    if ispyplot(Plotter)
        PyPlot = Plotter
        PyPlot.figure(1)
        PyPlot.clf()
        PyPlot.grid(true,color="gray", linestyle="dotted")
        PyPlot.title("Disk Limiting Current ($(disktype))")
        PyPlot.xlabel("sqrt(f/Hz)")
        PyPlot.ylabel("I/mA")
        frame=PyPlot.gca()
        #frame.axes.get_xaxis().set_ticks( [itick  for  itick=nfreq0:nfreq1])
        PyPlot.plot(sqrt.(freqs),res_idisk_sim./mA,label="FVM","ro-")
        PyPlot.plot(sqrt.(freqs),res_idisk_levich./mA,"g-",label="Levich")
        PyPlot.plot(sqrt.(freqs),res_idisk_newman./mA,"b-",label="Newman")
        PyPlot.legend(loc="upper left", fontsize="small", fancybox=false)
        PyPlot.tight_layout()
        PyPlot.savefig("levich-idisk.pdf")

        PyPlot.figure(2)
        PyPlot.clf()
        PyPlot.grid(true,color="gray", linestyle="dotted")
        PyPlot.title("Ring Limiting Current ($(disktype))")
        PyPlot.xlabel("sqrt(f/Hz)")
        PyPlot.ylabel("I/mA")
        frame=PyPlot.gca()
        #frame.axes.get_xaxis().set_ticks( [itick  for  itick=nfreq0:nfreq1])
        PyPlot.plot(sqrt.(freqs),res_iring_sim./mA,label="FVM", "ro-")
        PyPlot.plot(sqrt.(freqs),res_iring_levich./mA,"g-",label="Levich")
        PyPlot.legend(loc="upper left", fontsize="small", fancybox=false)
        PyPlot.tight_layout()
        PyPlot.savefig("levich-iring.pdf")

        PyPlot.figure(3)
        PyPlot.clf()
        PyPlot.grid(true,color="gray", linestyle="dotted")
        PyPlot.title("Collection Efficiency ($(disktype))")
        PyPlot.ylim(0.2,0.3)
        PyPlot.xlabel("sqrt(f/Hz)")
        PyPlot.ylabel("I/mA")
        frame=PyPlot.gca()
        #frame.axes.get_xaxis().set_ticks( [itick  for  itick=nfreq0:nfreq1])
        PyPlot.plot(sqrt.(freqs),res_coleff_albery,"g-",label="Albery")
        PyPlot.plot(sqrt.(freqs),res_coleff_sim,label="FVM", "ro-")
        PyPlot.legend(loc="upper right", fontsize="small", fancybox=false)
        PyPlot.tight_layout()
        PyPlot.savefig("levich-coleff.pdf")
    end
    return sum(res_idisk_sim)+sum(res_iring_sim)+sum(res_coleff_sim)

end
function runtests()
    result=main(maxfreq=4.0,verbose=false)
    @test result ≈ 0.5180855490694214 
end


function create_plots()
    main(plot_results=true,maxfreq=200.0)
end

end
