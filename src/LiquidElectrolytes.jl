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
VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public rexp, rlog"))

include("electrolyte.jl")
export ElectrolyteData, AbstractElectrolyteData
export dlcap0, chargedensity, chemical_potentials!, rrate, debyelength, chemical_potential, c0_barc
export showstruct, electrolyte, solventconcentration
export isincompressible, iselectroneutral
export chemical_potentials, electrochemical_potentials

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
export EquilibriumData, create_equilibrium_system, create_equilibrium_pp_system,  dlcapsweep_equi
VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public set_molarity!, apply_voltage!, update_derived!"))
VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public calc_cmol, calc_QBL, calc_φ, calc_p, ysum"))
include("equilibrium-supplement.jl")

include("pnpstokes.jl")
export PNPStokesSolver

include("RotatingElectrodes/RotatingElectrodes.jl")

end # module
