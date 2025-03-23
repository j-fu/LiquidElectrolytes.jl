"""
$(TYPEDEF)

Result data type for [`dlcapsweep`](@ref)

$(TYPEDFIELDS)

Access methods: [`voltages`](@ref), [`voltages_solutions`](@ref),  [`voltages_dlcaps`e](@ref)
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
        sys;
        electrolyte = electrolytedata(sys),
        inival = nothing,
        voltages = (-1:0.1:1) * ufac"V",
        molarity = nothing,
        δ = 1.0e-4,
        store_solutions = false,
        solver_kwargs...
    )
    (; ip, iϕ)= electrolyte
    if !isnothing(molarity)
        error("The molarity kwarg of dlcapsweep has been removed. Pass the molarity information with electrolyte.c_bulk.")
    end

    ranges = _splitz(voltages)

    electrolyte.ϕ_we = 0.0

    if isnothing(inival)
        inival = pnpunknowns(sys)
    end

    inival = solve(sys; inival, damp_initial = 0.1, solver_kwargs...)

    function show_error(u, δ)
        @show u[iϕ, 1:5]
        @show u[ip, 1:5]
        return @error "bailing out at δ=$(δ) ϕ_we=$(electrolyte.ϕ_we)V"
    end

    result_plus = DLCapSweepResult()
    result_minus = DLCapSweepResult()

    control = VoronoiFVM.SolverControl(; max_round = 3, tol_round = 1.0e-9,
                                       unorm = u -> wnorm(u, electrolyte.weights, Inf),
                                       rnorm = u -> wnorm(u, electrolyte.weights, 1),
                                       solver_kwargs...)
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
                electrolyte.ϕ_we = ϕ
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
                electrolyte.ϕ_we = ϕ + δ
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

