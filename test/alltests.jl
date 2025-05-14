using LiquidElectrolytes
using ExtendableFEM
using LinearAlgebra
using ExtendableGrids
using VoronoiFVM
using LessUnitful
using UUIDs
using ExampleJuggler
using ExplicitImports
using ExtendableFEM
using Markdown
using Aqua
using LiquidElectrolytes: aflux!, cflux!, sflux!

ExampleJuggler.verbose!(true)

@phconstants N_A
@unitfactors dm nm mol

thisproject = dirname(Base.active_project())

@testset "cdl0" begin
    ely = ElectrolyteData(c_bulk = fill(0.01 * mol / dm^3, 2) .|> unitfactor)
    @test dlcap0(ely) ≈ 0.22846691848825248
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
        acelldata = ElectrolyteData(; Γ_we = 1, Γ_bulk = 2, flux! = aflux!, κ, c_bulk = fill(molarity, 2))
        acell = PNPSystem(grid; bcondition, celldata = acelldata)
        aresult = dlcapsweep(acell; voltages, δ)
        avolts = aresult.voltages
        acaps = aresult.dlcaps

        μcelldata = ElectrolyteData(; Γ_we = 1, Γ_bulk = 2, flux! = sflux!, κ, c_bulk = fill(molarity, 2))
        μcell = PNPSystem(grid; bcondition, celldata = μcelldata)
        μresult = dlcapsweep(μcell; voltages, δ)
        μvolts = μresult.voltages
        μcaps = μresult.dlcaps

        ccelldata = ElectrolyteData(; Γ_we = 1, Γ_bulk = 2, flux! = cflux!, κ, c_bulk = fill(molarity, 2))
        ccell = PNPSystem(grid; bcondition, celldata = ccelldata)
        cresult = dlcapsweep(ccell; voltages, δ)
        cvolts = cresult.voltages
        ccaps = cresult.dlcaps


        function pbbcondition(f, u, bnode, data)
            (; Γ_we, Γ_bulk, ϕ_we, ip, iϕ) = data
            ## Dirichlet ϕ=ϕ_we at Γ_we
            boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)
            boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_bulk, value = data.ϕ_bulk)
            boundary_dirichlet!(f, u, bnode, species = ip, region = Γ_bulk, value = data.p_bulk)
        end

        pcelldata = ElectrolyteData(; Γ_we = 1, Γ_bulk = 2, κ, c_bulk = fill(molarity, 2))
        pcell = PBSystem(grid; bcondition = pbbcondition, celldata = pcelldata)
        presult = dlcapsweep(pcell; voltages, δ)
        pvolts = presult.voltages
        pcaps = presult.dlcaps

        @show norm((acaps - μcaps) ./ acaps)
        @show norm((acaps - ccaps) ./ acaps)
        @show norm((acaps - pcaps) ./ acaps)


        @test isapprox(dlcap0(acelldata), acaps[findfirst(x -> x ≈ 0, avolts)], rtol = 1.0e-2)
        @test isapprox(dlcap0(acelldata), ccaps[findfirst(x -> x ≈ 0, cvolts)], rtol = 1.0e-2)
        @test isapprox(dlcap0(acelldata), μcaps[findfirst(x -> x ≈ 0, μvolts)], rtol = 1.0e-2)
        @test isapprox(dlcap0(acelldata), pcaps[findfirst(x -> x ≈ 0, pvolts)], rtol = 1.0e-2)


        @test isapprox(acaps, μcaps, rtol = 1.0e-10)
        @test isapprox(acaps, ccaps, rtol = 1.0e-10)
        @test isapprox(acaps, pcaps, rtol = 1.0e-10)
    end
    xtest([0.0, 0.0])
    xtest([10.0, 10.0])

end


examples = [
    "Example101_DLCap.jl",
    "Example110_Fe23Cell.jl",
    #    "Example111_SurfaceKinetics.jl",
    "Example120_ORRCell.jl",
    "Example210_Leveque.jl",
    "Example211_Levich.jl",
    #    "Example220_SimpleCV.jl",
]

@testset "Examples" begin
    @testmodules(joinpath(@__DIR__, "..", "examples"), examples)
end


# Notebooks should "disable in file" the cell activating
# the docs enviroment
notebooks = [
    "EquilibriumCheck.jl",
    "ORR.jl",
    "PoissonBoltzmann.jl",
    "SurfaceKinetics_draft.jl",
    "BufferReactions.jl",
    "ElectroOsmosis.jl",
]


@testset "Notebooks" begin
    @testscripts(joinpath(@__DIR__, "..", "notebooks"), notebooks)
end


@testset "ExplicitImports" begin
    @test ExplicitImports.check_no_implicit_imports(LiquidElectrolytes, skip = (Base, Core, Markdown)) === nothing
    @test ExplicitImports.check_all_explicit_imports_are_public(LiquidElectrolytes) === nothing
    @test ExplicitImports.check_no_stale_explicit_imports(LiquidElectrolytes) === nothing
end

@testset "Aqua" begin
    Aqua.test_all(LiquidElectrolytes)
end

@testset "UndocumentedNames" begin
    if isdefined(Docs, :undocumented_names) # >=1.11
        @test isempty(Docs.undocumented_names(LiquidElectrolytes))
    end
end
