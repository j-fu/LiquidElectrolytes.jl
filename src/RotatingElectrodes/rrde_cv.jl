mutable struct CVContext
    reaction     # Disk/ring  reaction
    bulk_reaction # Reaction in bulk
    num_bulk_species  # number of bulk species
    set_voltage  # Set voltage - pre-caculate exponential terms
    rrde_data    # RRDE Data, see rrde_flux.jl  
    grid         # the grid
    rdcell       # the FVMSystem as rd cell 
    tfc          # the test functions
    nu           # Viscosity
    I_disk       # function calculating disk current from reacting species
    I_ring       # function calculating ring current from reacting species
    CVContext()=new()
end


"""
    Create default CV context for some reaction
"""
function CVContext(reaction,
                   disktype="E6R1";
                   set_voltage=()->null,
                   bulk_reaction=VoronoiFVM.nofunc,
                   nref=0,
                   num_species=1,
                   num_disk_species=0,
                   I_disk=(I)-> I[1],
                   I_ring=(I)-> I[1],
                   )
    cvcontext=CVContext()
    cvcontext.reaction=reaction
    cvcontext.bulk_reaction=bulk_reaction
    cvcontext.set_voltage=set_voltage
    cvcontext.rrde_data=RRDEData(ones(num_species))
    cvcontext.nu=1
    cvcontext.num_bulk_species=num_species
    cvcontext.I_disk=I_disk
    cvcontext.I_ring=I_ring
    geom=DiscSpec(disktype)
    cvcontext.grid=rrde_grid(geom;nref)


    function bstorage(f,u,node,data)
        if node.region==b_disk
            for ispec=num_species+1:num_species+num_disk_species
                f[ispec]=u[ispec]
            end
        end
    end

    function breaction(f,u,node,data)
        rrde_outflow(f,u,node,data)
        reaction(f,u,node,data)
    end


    
    if num_disk_species==0
        physics=VoronoiFVM.Physics(num_species=num_species+num_disk_species,
                                   data=cvcontext.rrde_data,
                                   reaction=bulk_reaction,
                                   flux=rrde_flux,
                                   breaction=breaction)
        cvcontext.rdcell=VoronoiFVM.DenseSystem(cvcontext.grid,physics)
    else
        physics=VoronoiFVM.Physics(num_species=num_species+num_disk_species,
                                   data=cvcontext.rrde_data,
                                   flux=rrde_flux,
                                   reaction=bulk_reaction,
                                   bstorage=bstorage,
                                   breaction=breaction)
        cvcontext.rdcell=VoronoiFVM.SparseSystem(cvcontext.grid,physics)
    end        
    
    for ispec=1:num_species
        enable_species!(cvcontext.rdcell,ispec,[1])
    end
    for ispec=num_species+1:num_species+num_disk_species
        enable_boundary_species!(cvcontext.rdcell,ispec,[b_disk])
    end

    
    
    factory=VoronoiFVM.TestFunctionFactory(cvcontext.rdcell)
    cvcontext.tfc=Vector{Any}(undef,b_max)
    cvcontext.tfc[b_disk]=testfunction(factory,[b_in,b_out,b_ring],[b_disk])
    cvcontext.tfc[b_ring]=testfunction(factory,[b_in,b_out,b_disk],[b_ring])
    cvcontext.tfc[b_in]=testfunction(factory,[b_out,b_disk,b_ring],[b_in])
    cvcontext.tfc[b_out]=testfunction(factory,[b_in,b_disk,b_ring],[b_out])

    
    cvcontext
end    

"""
    Set viscosity
"""
viscosity!(cvcontext,nu)= cvcontext.nu=nu

"""
    Set Diffusion coefficient
"""
diffusion!(cvcontext,ispec,D)= cvcontext.rrde_data.D[ispec]=D

"""
    Set Dirichlet boundary value
"""
VoronoiFVM.boundary_dirichlet!(cvcontext::CVContext,ispec,ibc,val)= VoronoiFVM.boundary_dirichlet!(cvcontext.rdcell,ispec,ibc,val)

"""
    Perform cv sweeps
"""
function run_cv(cvcontext;
                frequency=9,
                phi_min=0*V,
                phi_max=-1.6*V,
                scanrate=20*mV/s,
                num_periods=1,
                verbose=true,
                control=FVMNewtonControl()
                )
    @local_phconstants AvogadroConstant MolarGasConstant ElementaryCharge
    FaradayConstant = AvogadroConstant*ElementaryCharge
    @local_unitfactors mA
    
    
    period=2*abs(phi_max-phi_min)/scanrate
    sampling_times=[ i*period/2 for i=0:2*num_periods]

    omega=frequency*2.0*π
    @printf("---- rotation frequency: %.2f omega: %.2e  -------\n",frequency,omega)

    # adjustment of the grid relative to rotation frequency
    delta=1.61*(cvcontext.nu)^(1/6)*(cvcontext.rrde_data.D[1])^(1/3)/sqrt(omega)

    rescale_z!(cvcontext.grid,delta)
    update_grid!(cvcontext.rdcell)
    karman_data=KarmanData(omega,cvcontext.nu)
    
    fkarman(r,z)=fKarman(karman_data,r,z)
    update!(cvcontext.rrde_data,cvcontext.grid,fkarman)

    
    vdisk=[]
    iring=[]
    idisk=[]
    time_discretization=[]
    time=0.0

    function phi_cv(t)
        tper=t-period*floor(t/period)
        if tper<0.5*period
            phicv=phi_min+ 2*(tper/period)*(phi_max-phi_min)
        else
            phicv=phi_min+ 2*((period-tper)/period)*(phi_max-phi_min)
        end
        return phicv
    end

    
    solution=unknowns(cvcontext.rdcell)
    inival=unknowns(cvcontext.rdcell)
    inival.=0

    # solve for  inititial state 
    cvcontext.set_voltage(phi_cv(time))
    solution=solve(cvcontext.rdcell; inival)
    inival.=solution
    I_disk=VoronoiFVM.integrate(cvcontext.rdcell,cvcontext.tfc[b_disk],solution)*FaradayConstant
    I_ring=VoronoiFVM.integrate(cvcontext.rdcell,cvcontext.tfc[b_ring],solution)*FaradayConstant
    push!(vdisk,phi_cv(time))
    push!(iring,cvcontext.I_ring(I_ring))
    push!(idisk,cvcontext.I_disk(I_disk))
    push!(time_discretization,time)

    

        
    pmeter=ProgressUnknown("steps:")
    function pre(sol,time)
        cvcontext.set_voltage(phi_cv(time))
    end

    I_disk=zeros(cvcontext.num_bulk_species)
    I_disk_old=zeros(cvcontext.num_bulk_species)
    I_ring=[]
    di=0.0
    function delta(sys, solution, oldsolution, time,tstep)
        I_disk=VoronoiFVM.integrate(sys,cvcontext.tfc[b_disk],solution,oldsolution,tstep)*FaradayConstant
        I_ring=VoronoiFVM.integrate(sys,cvcontext.tfc[b_ring],solution,oldsolution,tstep)*FaradayConstant
        di=abs(cvcontext.I_disk(I_disk_old)-cvcontext.I_disk(I_disk))/mA
    end
    
        
    function post(solution, oldsolution, time,tstep)
        push!(vdisk,phi_cv(time))
        push!(iring,cvcontext.I_ring(I_ring))
        push!(idisk,cvcontext.I_disk(I_disk))
        push!(time_discretization,time)
        if verbose
            ProgressMeter.next!(pmeter,showvalues=[
                (:t,@sprintf("%.2f",time)),
                (:Δϕ,@sprintf("%.2f",phi_cv(time))),
                (:ΔI,@sprintf("%.2e",di)),
                (:Δt,@sprintf("%.2e",tstep)),
            ],valuecolor=:yellow)
        end
        I_disk_old=I_disk
    end

    solve(cvcontext.rdcell; inival, times=sampling_times, control=control, pre=pre,post=post,delta=delta)

    ProgressMeter.finish!(pmeter)            
    (frequency=frequency,vdisk=vdisk,iring=iring,idisk=idisk)
end


solvercontrol()=VoronoiFVM.NewtonControl()

CVResults()=(frequencies=[],irings=[],idisks=[],vdisks=[])

function add_cvresult!(cvresults,cvresult)
    push!(cvresults.frequencies,cvresult.frequency)
    push!(cvresults.vdisks,cvresult.vdisk)
    push!(cvresults.irings,cvresult.iring)
    push!(cvresults.idisks,cvresult.idisk)
end



plot_cvresults(cvresults;Plotter=nothing, file="cv.png",clear=true)=pyplot_cv(cvresults.frequencies,cvresults.vdisks,cvresults.idisks,cvresults.irings;
                                                                              Plotter, file=file,clear=clear)

