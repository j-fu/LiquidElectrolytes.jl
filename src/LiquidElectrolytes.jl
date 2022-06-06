module LiquidElectrolytes
using VoronoiFVM,ExtendableGrids
using DocStringExtensions
using Parameters
using ProgressLogging
using StaticArrays
using LinearAlgebra
using Unitful

function showstruct(io::IO,this)
    for name in fieldnames(typeof(this))
        println(io,"$(lpad(name,20)) = $(getfield(this,name))")
    end
end

include("units.jl")

@phconstants N_A e R ε_0 k_B
const F=N_A*e
const Mol=N_A
@siunits K  dm m s g nm Pa GPa V K L cm mA mol μF


export @siunits,@phconstants,@si_str


include("electrolyte.jl")
export ElectrolyteData, AbstractElectrolyteData
export Cdl0,chemical_potentials!, rrate,ldebye
export showstruct, rlog, electrolyte

include("pnpsystem.jl")
export PNPSystem
export pnpunknowns,electrolytedata,bulkbc


include("cells.jl")
export voltagesweep,doublelayercap


include("equilibrium.jl")
export EquilibriumData,L_debye, apply_voltage!,set_molarity!,update_derived!
export iφ,ip,iA,iC
export create_equilibrium_system, solve_equilibrium_system,create_equilibrium_pp_system
export calc_φ,calc_p, calc_cmol,calc_c0mol,calc_cnum, calc_QBL,ysum
export calc_Cdl,Cdl0, L_Debye
include("equilibrium-supplement.jl")


end # module
