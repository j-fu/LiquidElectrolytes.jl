module LiquidElectrolytes
using VoronoiFVM, ExtendableGrids
using DocStringExtensions
using ProgressLogging
using StaticArrays
using LinearAlgebra
using NLsolve
using Unitful, LessUnitful
import SciMLBase

using Base: @kwdef

function __init__()
    return LessUnitful.ensureSIBase()
end

include("tools.jl")
export RExp, RLog

include("electrolyte.jl")
export ElectrolyteData, AbstractElectrolyteData
export dlcap0, chargedensity, chemical_potentials!, rrate, debyelength, chemical_potential, c0_barc
export showstruct, electrolyte, solventconcentration
export isincompressible, iselectroneutral

include("pnpsystem.jl")
export PNPSystem
export pnpunknowns, electrolytedata, bulkbcondition

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public aflux, sflux, cflux"))

include("pbsystem.jl")
export PBSystem


include("cells.jl")
export ivsweep, dlcapsweep, currents, voltages_solutions, voltages_dlcaps, voltages_currents
export AbstractSimulationResult, DLCapSweepResult, IVSweepResult

include("equilibrium-pluto.jl")
export EquilibriumData, apply_voltage!, set_molarity!, update_derived!
#export iφ, ip, iA, iC
export create_equilibrium_system, solve_equilibrium_system, create_equilibrium_pp_system
export calc_φ, calc_p, calc_cmol, calc_c0mol, calc_cnum, calc_QBL, ysum, Cdl0
export dlcapsweep_equi
include("equilibrium-supplement.jl")

include("pnpstokes.jl")

include("RotatingElectrodes/RotatingElectrodes.jl")

end # module
