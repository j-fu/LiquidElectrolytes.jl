module LiquidElectrolytes
using VoronoiFVM, ExtendableGrids
using DocStringExtensions
using ProgressLogging
using StaticArrays
using LinearAlgebra
using NLsolve
using Unitful, LessUnitful

using SciMLBase
using ExtendableFEM

using Base: @kwdef

function __init__()
    LessUnitful.ensureSIBase()
end

myround(x; kwargs...) = round(x; kwargs...)
myround(s::Symbol; kwargs...) = s
myround(i::Int; kwargs...) = i
myround(b::Bool; kwargs...) = b

function showstruct(io::IO, this)
    for name in fieldnames(typeof(this))
        println(io, "$(lpad(name,20)) = $(myround.(getfield(this,name),sigdigits=5))")
    end
end

include("electrolyte.jl")
export ElectrolyteData, AbstractElectrolyteData
export dlcap0, chemical_potentials!, rrate, debyelength, chemical_potential, c0_barc
export showstruct, rlog, electrolyte, solventconcentration
export isincompressible, iselectroneutral

include("pnpsystem.jl")
export PNPSystem
export pnpunknowns, electrolytedata, bulkbcondition

include("pbsystem.jl")
export PBSystem


include("cells.jl")
export ivsweep, dlcapsweep, currents, voltages_solutions, voltages_dlcaps, voltages_currents
export AbstractSimulationResult, DLCapSweepResult, IVSweepResult

include("equilibrium-pluto.jl")
export EquilibriumData, apply_voltage!, set_molarity!, update_derived!
export iφ, ip, iA, iC
export create_equilibrium_system, solve_equilibrium_system, create_equilibrium_pp_system
export calc_φ, calc_p, calc_cmol, calc_c0mol, calc_cnum, calc_QBL, ysum, Cdl0
export dlcapsweep_equi
include("equilibrium-supplement.jl")


include("flowsolver.jl")
include("pnpstokes.jl")

end # module
