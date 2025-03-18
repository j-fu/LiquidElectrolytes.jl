module RotatingElectrodes
using DocStringExtensions: DocStringExtensions, SIGNATURES
using ExtendableGrids: ExtendableGrids, YCoordinates, bfacemask!, geomspace,
glue, linspace
using LessUnitful: LessUnitful, @local_phconstants, @local_unitfactors, @ph_str
using Printf:  @printf, @sprintf
using ProgressMeter: ProgressMeter, ProgressUnknown
using VoronoiFVM: VoronoiFVM, bfacevelocities, circular_symmetric!, coordinates,
edgevelocities, enable_boundary_species!, enable_species!,
    solve, testfunction, unknowns, update_grid!


include("asymptotics.jl")
export idisk_levich, idisk_newman, iring_levich, coleff_albery, dlayer

include("discspec.jl")
export DiscSpec

include("rrde_grid.jl")
export rrde_grid, rescale_z!
export b_inert, b_in, b_out, b_disk, b_ring, b_symm, b_max


include("karman.jl")
export KarmanData, KarmanData!, fKarman

include("rrde_flux.jl")
export RRDEData, rrde_flux, update!, rrde_outflow

include("rrde_cv.jl")
export CVContext, CVResults
export viscosity!, diffusion!, solvercontrol
export run_cv, add_cvresult!, plot_cvresults

include("pyplot.jl")
export pyplot_cv

end
