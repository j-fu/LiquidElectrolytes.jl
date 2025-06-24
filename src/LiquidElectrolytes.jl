"""
$(README)
"""
module LiquidElectrolytes
using Base: @kwdef
using DocStringExtensions: DocStringExtensions, TYPEDEF, TYPEDFIELDS, README
using ExtendableGrids: ExtendableGrids, ExtendableGrid, num_nodes, num_cellregions
using LessUnitful: LessUnitful, @local_phconstants, @local_unitfactors, @ph_str, @ufac_str
using InteractiveUtils: InteractiveUtils
using Markdown: @md_str
using Printf: Printf, @printf, @sprintf
using ProgressMeter: ProgressMeter
using ProgressLogging: @withprogress, @logprogress
using SciMLBase: SciMLBase, solve!
using VoronoiFVM: VoronoiFVM, TransientSolution, enable_boundary_species!, enable_species!, solve, testfunction, unknowns
using VoronoiFVM: boundary_dirichlet!, fbernoulli_pm, SolverControl
import VoronoiFVM
using LinearAlgebra: LinearAlgebra
using PreallocationTools: DiffCache, get_tmp

function __init__()
    return LessUnitful.ensureSIBase()
end

include("utils.jl")
export RExp, RLog

include("electrolyte.jl")
export AbstractElectrochemicalSystem, ElectrolyteData, AbstractElectrolyteData, update_derived!
export dlcap0, chargedensity, chemical_potentials!, rrate, debyelength, chemical_potential, c0_barc, solventconcentration
export isincompressible, iselectroneutral
export chemical_potentials, electrochemical_potentials

include("celldata.jl")

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse(
"""
public AbstractCellData, electrolytes, working_electrode, bulk_electrode,
   norm_weights, working_electrode_voltage, working_electrode_voltage!, 
   pressure_index, voltage_index, check_celldata
"""))

include("pnpsystem.jl")
export PNPSystem
export electrolytedata, celldata, bulkbcondition, unknowns

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public act_flux!, μex_flux!, cent_flux!, DGML_gamma!"))

include("pbsystem.jl")
export PBSystem


include("results.jl")
export AbstractSimulationResult, voltages, currents, voltages_solutions, voltages_dlcaps, voltages_currents
include("dlcapsweep.jl")
export dlcapsweep, DLCapSweepResult
include("ivsweep.jl")
export ivsweep, IVSweepResult
include("cvsweep.jl")
export SawTooth, cvsweep, CVSweepResult, period

include("pnpstokes.jl")
export PNPStokesSolver

include("RotatingElectrodes/RotatingElectrodes.jl")

end # module
