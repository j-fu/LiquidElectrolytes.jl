"""
    LiquidElectrolytes

$(read(joinpath(@__DIR__, "..", "README.md"), String))
"""
module LiquidElectrolytes
using Base: @kwdef
using DocStringExtensions: DocStringExtensions,  TYPEDEF, TYPEDFIELDS
using ExtendableGrids: ExtendableGrids, num_nodes
using LessUnitful: LessUnitful, @local_phconstants, @local_unitfactors, @ph_str, @ufac_str
using InteractiveUtils: InteractiveUtils
using Markdown: @md_str
using Printf: Printf, @printf, @sprintf
using ProgressMeter: ProgressMeter
using ProgressLogging: @withprogress, @logprogress
using SciMLBase: SciMLBase, solve!
using VoronoiFVM: VoronoiFVM, TransientSolution, enable_boundary_species!, solve, testfunction, unknowns
using VoronoiFVM: boundary_dirichlet!, fbernoulli_pm, SolverControl
using LinearAlgebra: LinearAlgebra

function __init__()
    return LessUnitful.ensureSIBase()
end

include("utils.jl")
export RExp, RLog
VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public rexp, rlog"))

include("electrolyte.jl")
export ElectrolyteData, AbstractElectrolyteData
export dlcap0, chargedensity, chemical_potentials!, rrate, debyelength, chemical_potential, c0_barc, solventconcentration
export isincompressible, iselectroneutral
export chemical_potentials, electrochemical_potentials

include("pnpsystem.jl")
export PNPSystem
export pnpunknowns, electrolytedata, bulkbcondition

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public aflux, sflux, cflux"))

include("pbsystem.jl")
export PBSystem


include("results.jl")
export AbstractSimulationResult, voltages, currents, voltages_solutions, voltages_dlcaps, voltages_currents
include("dlcapsweep.jl")
export dlcapsweep,  DLCapSweepResult
include("ivsweep.jl")
export ivsweep, IVSweepResult
include("cvsweep.jl")
export SawTooth, cvsweep, CVSweepResult, period

include("pnpstokes.jl")
export PNPStokesSolver

include("RotatingElectrodes/RotatingElectrodes.jl")

end # module
