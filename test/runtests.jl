using Test
using LiquidElectrolytes
using LinearAlgebra
using ExtendableGrids
using VoronoiFVM
using LessUnitful
using Pluto

@phconstants N_A
@unitfactors dm nm mol

@testset "cdl0" begin
    ely=ElectrolyteData(c_bulk=fill(0.01*mol/dm^3,2).|>unitfactor)
    @test dlcap0(ely)≈ 0.22846691848825248
    edata=EquilibriumData()
    LiquidElectrolytes.set_molarity!(edata,0.01)
    edata.χ=78.49-1
    @test dlcap0(edata)|>unitfactor≈ 0.22846691848825248
    @test dlcap0(EquilibriumData(ely))≈ 0.22846691848825248
end


@testset "dlcap" begin
    
    
    function bcondition(f,u,bnode,data::ElectrolyteData)
        (; iϕ,Γ_we,ϕ_we) = data
        boundary_dirichlet!(f,u,bnode,species=iϕ,region=Γ_we,value=ϕ_we)
        bulkbcondition(f,u,bnode,data)
    end


    hmin=1.0e-1*nm
    hmax=1.0*nm
    L=20.0*nm
    X=geomspace(0,L,hmin,hmax)
    molarity=0.1*ufac"M"
    voltages=collect(-1:0.01:1)*ufac"V"
    δ=1.0e-4
    
    grid=simplexgrid(X)
    function xtest(κ)
        acelldata=ElectrolyteData(;Γ_we=1, Γ_bulk=2,  scheme=:act,κ)
        acell=PNPSystem(grid;bcondition,celldata=acelldata)
        aresult=dlcapsweep(acell;voltages,molarity,δ)
        avolts=aresult.voltages
        acaps=aresult.dlcaps
        
        μcelldata=ElectrolyteData(;Γ_we=1, Γ_bulk=2, scheme=:μex,κ)
        μcell=PNPSystem(grid;bcondition,celldata=μcelldata)
        μresult=dlcapsweep(μcell;voltages,molarity,δ)
        μvolts=μresult.voltages
        μcaps=μresult.dlcaps
        
        ccelldata=ElectrolyteData(;Γ_we=1, Γ_bulk=2, scheme=:cent,κ)
        ccell=PNPSystem(grid;bcondition,celldata=ccelldata)
        cresult=dlcapsweep(ccell;voltages,molarity,δ)
        cvolts=cresult.voltages
        ccaps=cresult.dlcaps
        
        
        function pbbcondition(f,u,bnode,data)
	    (;Γ_we,Γ_bulk,ϕ_we) = data
       	    iϕ,ip=1,2 
	    ## Dirichlet ϕ=ϕ_we at Γ_we
	    boundary_dirichlet!(f,u,bnode,species=iϕ,region=Γ_we,value=ϕ_we)
	    boundary_dirichlet!(f,u,bnode,species=iϕ,region=Γ_bulk,value=data.ϕ_bulk)
	    boundary_dirichlet!(f,u,bnode,species=ip,region=Γ_bulk,value=data.p_bulk)
        end
        
        pcelldata=ElectrolyteData(;Γ_we=1, Γ_bulk=2, scheme=:cent,κ)
        pcell=PBSystem(grid;bcondition=pbbcondition,celldata=pcelldata)
        presult=dlcapsweep(pcell;inival=unknowns(pcell),voltages,molarity,δ, iϕ=1)
        pvolts=presult.voltages
        pcaps=presult.dlcaps
        
        ecell=create_equilibrium_system(grid,EquilibriumData(acelldata))
        evolts,ecaps=dlcapsweep_equi(ecell,vmax=1.0,molarity=0.1,δV=1.0e-4,nsteps=101)
        
        @show norm((acaps-ecaps)./ecaps)
        @show norm((acaps-μcaps)./acaps)
        @show norm((acaps-ccaps)./acaps)
        @show norm((acaps-pcaps)./acaps)
        
        # @show dlcap0(acelldata)
        # @show acaps[findfirst(x->x≈0,avolts)]
        # @show ecaps[findfirst(x->x≈0,evolts)]
        # @show ccaps[findfirst(x->x≈0,cvolts)]
        # @show μcaps[findfirst(x->x≈0,μvolts)]
        # @show pcaps[findfirst(x->x≈0,pvolts)]
        
        @test dlcap0(acelldata) ≈ dlcap0(EquilibriumData(acelldata))
        
        @test isapprox(dlcap0(acelldata),acaps[findfirst(x->x≈0,avolts)], rtol=1.0e-2)
        @test isapprox(dlcap0(acelldata),ecaps[findfirst(x->x≈0,evolts)], rtol=1.0e-2)
        @test isapprox(dlcap0(acelldata),ccaps[findfirst(x->x≈0,cvolts)], rtol=1.0e-2) 
        @test isapprox(dlcap0(acelldata),μcaps[findfirst(x->x≈0,μvolts)], rtol=1.0e-2)
        @test isapprox(dlcap0(acelldata),pcaps[findfirst(x->x≈0,pvolts)], rtol=1.0e-2)
        

        @test isapprox(acaps,ecaps,rtol=1.0e-3)
        @test isapprox(acaps,μcaps,rtol=1.0e-10)
        @test isapprox(acaps,ccaps,rtol=1.0e-10)
        @test isapprox(acaps,pcaps,rtol=1.0e-10)
    end
    xtest([0.0,0.0])
    xtest([10.0,10.0])
    
end

function run_notebook_in_current_environment(notebookname)
    # Prevent Pluto from calling into registry update
    Pluto.PkgCompat._updated_registries_compat[]=true

    session = Pluto.ServerSession();
    session.options.server.disable_writing_notebook_files=true
    session.options.server.show_file_system=false
    session.options.server.launch_browser=false
    session.options.server.dismiss_update_notification=true
    session.options.evaluation.capture_stdout=false
    session.options.evaluation.workspace_use_distributed=false # this makes it fast

    wd=pwd()
    t= @elapsed notebook = Pluto.SessionActions.open(session, notebookname; run_async=false)
    @info "notebook executed in $(round(t,sigdigits=4)) seconds"
    cd(wd)
    errored=false
    for c in notebook.cells
        if c.errored
            errored=true
            @error "Error in  $(c.cell_id): $(c.output.body[:msg])\n $(c.code)"
        end
    end
    !errored
end


notebooks=["ORR.jl"]

if VERSION>v"1.8" && !(VERSION>v"1.9.99")
    ENV["PLUTO_DEVELOP"]=true
    ENV["PLUTO_CI"]=true
    @testset "notebooks" begin
        for notebook in notebooks
            @info "notebook $(notebook):"
            @test run_notebook_in_current_environment(joinpath(@__DIR__,"..","notebooks",notebook))
            @info "notebook $(notebook) ok"
        end
    end
end
