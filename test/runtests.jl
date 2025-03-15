using Test
using LiquidElectrolytes
using ExtendableFEM
using LinearAlgebra
using ExtendableGrids
using VoronoiFVM
using LessUnitful
using UUIDs
using ExampleJuggler
ExampleJuggler.verbose!(true)

@phconstants N_A
@unitfactors dm nm mol

thisproject = dirname(Base.active_project())

@testset "cdl0" begin
    ely = ElectrolyteData(c_bulk = fill(0.01 * mol / dm^3, 2) .|> unitfactor)
    @test dlcap0(ely) ≈ 0.22846691848825248
    edata = EquilibriumData()
    LiquidElectrolytes.set_molarity!(edata, 0.01)
    edata.χ = 78.49 - 1
    @test dlcap0(edata) |> unitfactor ≈ 0.22846691848825248
    @test dlcap0(EquilibriumData(ely)) ≈ 0.22846691848825248
end


@testset "dlcap" begin


    function bcondition(f, u, bnode, data::ElectrolyteData)
        (; iϕ, Γ_we, ϕ_we) = data
        boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)
        bulkbcondition(f, u, bnode, data)
    end


    hmin = 1.0e-1 * nm
    hmax = 1.0 * nm
    L = 20.0 * nm
    X = geomspace(0, L, hmin, hmax)
    molarity = 0.1 * ufac"M"
    voltages = collect(-1:0.005:1) * ufac"V"
    δ = 1.0e-4

    grid = simplexgrid(X)
    function xtest(κ)
        acelldata = ElectrolyteData(; Γ_we = 1, Γ_bulk = 2, scheme = :act, κ, c_bulk = fill(molarity , 2))
        acell = PNPSystem(grid; bcondition, celldata = acelldata)
        aresult = dlcapsweep(acell; voltages, δ)
        avolts = aresult.voltages
        acaps = aresult.dlcaps

        μcelldata = ElectrolyteData(; Γ_we = 1, Γ_bulk = 2, scheme = :μex, κ, c_bulk = fill(molarity , 2))
        μcell = PNPSystem(grid; bcondition, celldata = μcelldata)
        μresult = dlcapsweep(μcell; voltages, δ)
        μvolts = μresult.voltages
        μcaps = μresult.dlcaps

        ccelldata = ElectrolyteData(; Γ_we = 1, Γ_bulk = 2, scheme = :cent, κ, c_bulk = fill(molarity, 2))
        ccell = PNPSystem(grid; bcondition, celldata = ccelldata)
        cresult = dlcapsweep(ccell; voltages, δ)
        cvolts = cresult.voltages
        ccaps = cresult.dlcaps


        function pbbcondition(f, u, bnode, data)
            (; Γ_we, Γ_bulk, ϕ_we) = data
            iϕ, ip = 1, 2
            ## Dirichlet ϕ=ϕ_we at Γ_we
            boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)
            boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_bulk, value = data.ϕ_bulk)
            boundary_dirichlet!(f, u, bnode, species = ip, region = Γ_bulk, value = data.p_bulk)
        end

        pcelldata = ElectrolyteData(; Γ_we = 1, Γ_bulk = 2, scheme = :cent, κ, c_bulk = fill(molarity , 2))
        pcell = PBSystem(grid; bcondition = pbbcondition, celldata = pcelldata)
        presult = dlcapsweep(pcell; inival = unknowns(pcell, inival = 0), voltages, δ, iϕ = 1)
        pvolts = presult.voltages
        pcaps = presult.dlcaps

        ecelldata=EquilibriumData(acelldata)
        ecell = create_equilibrium_system(grid,ecelldata)
        LiquidElectrolytes.set_molarity!(ecelldata, 0.1)
        evolts, ecaps = dlcapsweep_equi(ecell, vmax = 1.0, δV = 1.0e-4, nsteps = 201)

        @show norm((acaps - ecaps) ./ ecaps)
        @show norm((acaps - μcaps) ./ acaps)
        @show norm((acaps - ccaps) ./ acaps)
        @show norm((acaps - pcaps) ./ acaps)

        # @show dlcap0(acelldata)
        # @show acaps[findfirst(x->x≈0,avolts)]
        # @show ecaps[findfirst(x->x≈0,evolts)]
        # @show ccaps[findfirst(x->x≈0,cvolts)]
        # @show μcaps[findfirst(x->x≈0,μvolts)]
        # @show pcaps[findfirst(x->x≈0,pvolts)]

        @test dlcap0(acelldata) ≈ dlcap0(EquilibriumData(acelldata))

        @test isapprox(dlcap0(acelldata), acaps[findfirst(x -> x ≈ 0, avolts)], rtol = 1.0e-2)
        @test isapprox(dlcap0(acelldata), ecaps[findfirst(x -> x ≈ 0, evolts)], rtol = 1.0e-2)
        @test isapprox(dlcap0(acelldata), ccaps[findfirst(x -> x ≈ 0, cvolts)], rtol = 1.0e-2)
        @test isapprox(dlcap0(acelldata), μcaps[findfirst(x -> x ≈ 0, μvolts)], rtol = 1.0e-2)
        @test isapprox(dlcap0(acelldata), pcaps[findfirst(x -> x ≈ 0, pvolts)], rtol = 1.0e-2)


        @test isapprox(acaps, ecaps, rtol = 1.0e-3)
        @test isapprox(acaps, μcaps, rtol = 1.0e-10)
        @test isapprox(acaps, ccaps, rtol = 1.0e-10)
        @test isapprox(acaps, pcaps, rtol = 1.0e-10)
    end
    xtest([0.0, 0.0])
    xtest([10.0, 10.0])

end

notebooks = [
    "DLCap.jl",
    "ORR.jl",
    "Equilibrium1D.jl",
    "PoissonBoltzmann.jl",
    "SurfaceKinetics_draft.jl",
    "BufferReactions.jl",
    "ElectroOsmosis.jl",
]

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"))
Pkg.develop(path = joinpath(@__DIR__, ".."))
Pkg.instantiate()
using LiquidElectrolytes
using ExtendableFEM

@testset "Notebooks" begin
    @testscripts(joinpath(@__DIR__, "..", "notebooks"), notebooks)
end
