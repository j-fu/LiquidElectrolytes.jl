"""
$(TYPEDEF)

Result data type for [`dlcapsweep`](@ref)

$(TYPEDFIELDS)

Access methods: [`voltages`](@ref), [`voltages_solutions`](@ref),  [`voltages_dlcaps`](@ref)
"""
Base.@kwdef struct DLCapSweepResult <: AbstractSimulationResult
    """
    Vector of voltages
    """
    voltages::Vector = zeros(0)

    """
    Vector of double layer capacitances
    """
    dlcaps::Vector = zeros(0)

    """
    Vector of solutions
    """
    solutions::Vector = []
end


"""
    function dlcapsweep(
            sys;
            electrolyte = electrolytedata(sys),
            inival = nothing,
            voltages = (-1:0.1:1) * ufac"V",
            δ = 1.0e-4,
            store_solutions = false,
            solver_kwargs...
        )
Calculate double layer capacitances for voltages given in `voltages`.
Returns a [`DLCapSweepResult`](@ref).

Assumptions:
- Only one double layer in the system - close to working electrode.
"""
function dlcapsweep(
        esys::AbstractElectrochemicalSystem;
        cdata = celldata(esys),
        inival = nothing,
        voltages = (-1:0.1:1) * ufac"V",
        molarity = nothing,
        δ = 1.0e-4,
        store_solutions = false,
        solver_kwargs...
    )
    update_derived!(cdata)
    iϕ = voltage_index(cdata)
    sys = esys.vfvmsys
    if !isnothing(molarity)
        error("The molarity kwarg of dlcapsweep has been removed. Pass the molarity information with electrolyte.c_bulk.")
    end

    ranges = splitz(voltages)

    working_electrode_voltage!(cdata, 0.0)

    if isnothing(inival)
        inival = unknowns(esys)
    end

    inival = solve(sys; inival, damp_initial = 0.1, solver_kwargs...)

    function show_error(u, δ)
        return @error "bailing out at δ=$(δ) ϕ_we=$(working_electrode_voltage(cdata))V"
    end

    result_plus = DLCapSweepResult()
    result_minus = DLCapSweepResult()

    control = VoronoiFVM.SolverControl(;
        max_round = 3, tol_round = 1.0e-9,
        unorm = u -> wnorm(u, norm_weights(cdata), Inf),
        rnorm = u -> wnorm(u, norm_weights(cdata), 1),
        solver_kwargs...
    )
    for range in ranges
        sol = inival
        success = true
        allprogress = length(range)
        if range[end] > range[1]
            result = result_plus
        else
            result = result_minus
        end
        ϕprogress = 0
        @withprogress name = "sweep $(range[1])V -> $(range[end])V" for ϕ in range
            try
                working_electrode_voltage!(cdata, ϕ)
                sol = solve(sys; inival = sol, control)
            catch e
                println(e)
                show_error(sol, 0)
                success = false
            end


            if !success
                break
            end
            store_solutions ? push!(result.solutions, sol) : nothing
            Q = VoronoiFVM.integrate(sys, sys.physics.reaction, sol)

            try
                working_electrode_voltage!(cdata, ϕ + δ)
                sol = solve(sys; inival = sol, control)
            catch e
                println(e)
                show_error(sol, δ)
                success = false
            end
            if !success
                break
            end
            Qδ = VoronoiFVM.integrate(sys, sys.physics.reaction, sol)

            cdl = (Qδ[iϕ] - Q[iϕ]) / δ

            push!(result.voltages, ϕ)
            push!(result.dlcaps, cdl)
            ϕprogress += 1
            @logprogress ϕprogress / allprogress V = ϕprogress
        end
    end

    return DLCapSweepResult(
        voltages = vcat(reverse(result_minus.voltages), result_plus.voltages),
        dlcaps = vcat(reverse(result_minus.dlcaps), result_plus.dlcaps),
        solutions = vcat(reverse(result_minus.solutions), result_plus.solutions)
    )
end
