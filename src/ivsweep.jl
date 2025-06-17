"""
    $(TYPEDEF)

Result data type for [`ivsweep`](@ref)

$(TYPEDFIELDS)

Access methods: [`voltages`](@ref), [`currents`](@ref),  [`voltages_currents`](@ref), [`voltages_solutions`](@ref)
"""
Base.@kwdef mutable struct IVSweepResult <: AbstractSimulationResult
    """
    Vector of voltages
    """
    voltages = zeros(0)

    """
    Working electrode molar reaction rates
    """
    j_reaction = []

    """
    Working electrode molar fluxes
    """
    j_we = []

    """
    Bulk molar fluxes
    """
    j_bulk = []

    """
    Vector of solutions
    """
    solutions = []
end

"""
     ivsweep(
          sys;
          voltages = (-0.5:0.1:0.5) * ufac"V",
          store_solutions = false,
          solver_kwargs...,
          )

Steady state voltammetry.

Calculate molar reaction rates and bulk flux rates for each voltage in `voltages`.
Returns a [`IVSweResult`](@ref).
"""
function ivsweep(
        esys::AbstractElectrochemicalSystem;
        cdata = celldata(esys),
        voltages = (-0.5:0.1:0.5) * ufac"V",
        store_solutions = false,
        solver_kwargs...
    )
    update_derived!(cdata)
    sys = esys.vfvmsys
    ranges = splitz(voltages)
    F = ph"N_A" * ph"e"
    factory = VoronoiFVM.TestFunctionFactory(sys)
    tf_bulk = testfunction(factory, [working_electrode(cdata)], [bulk_electrode(cdata)])
    tf_we = testfunction(factory, [bulk_electrode(cdata)], [working_electrode(cdata)])

    iplus = []
    iminus = []
    fplus = []
    fminus = []
    vminus = zeros(0)
    vplus = zeros(0)
    sminus = []
    splus = []

    working_electrode_voltage!(cdata, 0)
    control = SolverControl(;
        verbose = true,
        handle_exceptions = true,
        Δp_min = 1.0e-3,
        Δp = 1.0e-2,
        Δp_grow = 1.2,
        #          Δu_opt = 1.0e-2,
        unorm = u -> wnorm(u, norm_weights(cdata), Inf),
        rnorm = u -> wnorm(u, norm_weights(cdata), 1),
        solver_kwargs...
    )

    @info "Solving for 0V..."
    inival = solve(sys; inival = unknowns(esys), control)

    result_plus = IVSweepResult()
    result_minus = IVSweepResult()

    allprogress = voltages[end] - voltages[begin]
    ϕprogress = 0
    for range in ranges
        @info "IV sweep from $(range[1])V to $(range[end])V..."
        dir = range[end] > range[1] ? 1 : -1

        if dir > 0
            result = result_plus
        else
            result = result_minus
        end
        psol = nothing
        @withprogress begin
            function pre(sol, ϕ)
                working_electrode_voltage!(cdata,dir * ϕ)
            end

            function post(sol, oldsol, ϕ, Δϕ)
                I = -VoronoiFVM.integrate(sys, sys.physics.breaction, sol; boundary = true)
                I_react = I[:, working_electrode(cdata)]
                I_bulk = -VoronoiFVM.integrate(sys, tf_bulk, sol)
                I_we = -VoronoiFVM.integrate(sys, tf_we, sol)
                push!(result.voltages, working_electrode_voltage(cdata))
                push!(result.j_we, I_we)
                push!(result.j_bulk, I_bulk)
                push!(result.j_reaction, I_react)
                ϕprogress += abs(Δϕ)
                @logprogress ϕprogress / allprogress
            end

            function delta(sys, u, v, t, Δt)
                n = wnorm(u - v, norm_weights(cdata), Inf)
            end

            psol = solve(
                sys;
                inival,
                embed = dir * range,
                control,
                pre,
                post,
                delta,
                store_all = store_solutions
            )
        end
        if store_solutions
            result.solutions = psol.u
        end
        if dir == -1
            popfirst!(result.j_we)
            popfirst!(result.j_reaction)
            popfirst!(result.j_bulk)
            popfirst!(result.voltages)
            if store_solutions
                popfirst!(result.solutions)
            end
        end
    end

    return IVSweepResult(;
        voltages = vcat(reverse(result_minus.voltages), result_plus.voltages),
        j_reaction = vcat(reverse(result_minus.j_reaction), result_plus.j_reaction),
        j_bulk = vcat(reverse(result_minus.j_bulk), result_plus.j_bulk),
        j_we = vcat(reverse(result_minus.j_we), result_plus.j_we),
        solutions = vcat(reverse(result_minus.solutions), result_plus.solutions)
    )
end
