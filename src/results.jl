"""
    $(TYPEDEF)

Abstract supertype for simulation results.
"""
abstract type AbstractSimulationResult end


"""
    voltages_solutions(result)

Return a [`TransientSolution`](https://j-fu.github.io/VoronoiFVM.jl/stable/solutions/#Transient-solution) `tsol`
containing voltages (in `tsol.t`) and the corresponding stationary solutions (in `tsol.u`).
"""
voltages_solutions(r::AbstractSimulationResult) = TransientSolution(r.solutions, r.voltages)

"""
    voltages(result)

Return vector of voltages of simulation result.
"""
voltages(r::AbstractSimulationResult) = r.voltages


"""
    voltages_currents(result,ispec; electrode)

Voltage- working electrode current curve for species as [`DiffEqArray`](https://docs.sciml.ai/RecursiveArrayTools/stable/array_types/#RecursiveArrayTools.DiffEqArray)
"""
function voltages_currents(r::AbstractSimulationResult, ispec; electrode = :we)
    RecursiveArrayTools.DiffEqArray(currents(r, ispec; electrode), r.voltages)
end


"""
    currents(result,ispec; electrode)

Vector of electrode currents for species `ispec`.
Electrode can be: 
- `:bulk`: calculate current at bulk boundary conditions
- `:we`: calculate current electrode interface
- `:reaction`: calculate current from reaction term
"""
function currents(r::AbstractSimulationResult, ispec; electrode = :we)
    F = ph"N_A" * ph"e"
    return if electrode == :we
        [F * j[ispec] for j in r.j_we]
    elseif electrode == :bulk
        [F * j[ispec] for j in r.j_bulk]
    elseif electrode == :reaction
        return [F * j[ispec] for j in r.j_reaction]
    else
        error("no such electrode")
    end
end

"""
    voltages_dlcaps(result)

Double layer capacitance curve as [`DiffEqArray`](https://docs.sciml.ai/RecursiveArrayTools/stable/array_types/#RecursiveArrayTools.DiffEqArray)
"""
voltages_dlcaps(r::AbstractSimulationResult) = RecursiveArrayTools.DiffEqArray(r.dlcaps, r.voltages)
