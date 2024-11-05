module RotatingElectrodes
using Printf


using ExtendableGrids
using VoronoiFVM
using ProgressMeter
using PyPlot
using DocStringExtensions

include("units.jl")
using .UnitsConstants



include("asymptotics.jl")
export idisk_levich,idisk_newman,iring_levich,coleff_albery,dlayer

include("discspec.jl")
export DiscSpec

include("rrde_grid.jl")
export rrde_grid_zequi,rescale_z!
export b_inert,b_in,b_out,b_disk,b_ring,b_symm,b_max


include("karman.jl")
export KarmanData,KarmanData!,fKarman

include("rrde_flux.jl")
export RRDEData,rrde_flux,update!,rrde_outflow

include("rrde_cv.jl")
export CVContext,CVResults
export viscosity!,diffusion!,solvercontrol
export run_cv,add_cvresult!,plot_cvresults

include("pyplot.jl")
export pyplot_cv

end

