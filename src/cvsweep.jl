"""
    $(TYPEDEF)

Callable (functor) struct defining sawtooth function for cyclic voltammetry.

`(::SawTooth)(t)` returns the voltage to be applied at moment t.

$(TYPEDFIELDS)
"""
Base.@kwdef struct SawTooth
    "Minimum voltage. Default: -1V."
    vmin::Float64 = -1 * ufac"V"

    "Maximum voltage. Default: 1V."
    vmax::Float64 = 1 * ufac"V"

    "Scan rate. Default: 0.1V/s."
    scanrate::Float64 = 0.1 * ufac"V/s"

    "Start voltage. Default: 0V"
    vstart::Float64 = 0 * ufac"V"

    "Start time. Default: 0s"
    tstart::Float64 = 0 * ufac"s"

    "Initial scan direction. Default: true"
    scanup::Bool = true
end

function (this::SawTooth)(t::Number)
    (; vmin, vmax, scanrate, vstart, tstart, scanup) = this
    fullperiod = period(this)
    halfperiod = fullperiod / 2
    if scanup
        t = t - tstart + (vstart - vmin) / scanrate
    else
        t = t - tstart + (vmax - vstart) / scanrate + halfperiod
    end
    iperiod = Int(floor(t / fullperiod))
    t = t - iperiod * fullperiod
    if t < halfperiod
        v = vmin + t * scanrate
    else
        t = t - halfperiod
        v = vmax - t * scanrate
    end
    return v
end

"""
    period(st::SawTooth)

Return length of CV period 
"""
function period(st::SawTooth)
    (; vmin, vmax, scanrate, vstart, tstart, scanup) = st
    return 2 * (vmax - vmin) / scanrate
end

"""
    $(TYPEDEF)

Result type for [`cvsweep`](@ref).

Fields:
$(TYPEDFIELDS)

Access methods: [`voltages`](@ref), [`currents`](@ref)  [`voltages_currents`](@ref)

"""
Base.@kwdef mutable struct CVSweepResult <: AbstractSimulationResult
    """
    Vector of voltages
    """
    voltages = zeros(0)

    """
    Vector of times
    """
    times = zeros(0)

    """
    Working electrode molar reaction rates
    """
    j_reaction = []

    """
    Bulk molar fluxes
    """
    j_bulk = []

    """
    Working electrode molar fluxes
    """
    j_we = []

    """
    Transient solution
    """
    tsol = nothing
end

"""
     cvsweep(
          sys;
          voltages=SawTooth(),
          periods=1,
          solver_kwargs...,
          )

Cyclic voltammetry sweep. By default, `voltages` given as a [`SawTooth`](@ref) function.

Calculate molar reaction rates and working electrode flux rates.
Returns a [`CVSweepResult`](@ref).
"""
function cvsweep(
        esys::AbstractElectrochemicalSystem;
        cdata = celldata(esys),
        voltages = SawTooth(),
        nperiods = 1,
        store_solutions = false,
        solver_kwargs...
    )
    update_derived!(cdata)
    sys = esys.vfvmsys
    F = ph"N_A" * ph"e"
    factory = VoronoiFVM.TestFunctionFactory(sys)
    tf_we = testfunction(factory, [bulk_electrode(cdata)], [working_electrode(cdata)])
    tf_bulk = testfunction(factory, [working_electrode(cdata)], [bulk_electrode(cdata)])

    working_electrode_voltage!(cdata, voltages(0))
    control = SolverControl(;
        verbose = "",
        handle_exceptions = true,
        damp_initial = 1,
        Δu_opt = 0.05,
        Δt_min = 1.0e-4 * period(voltages),
        Δt_max = 1.0e-2 * period(voltages),
        Δt = 1.0e-3 * period(voltages),
        Δt_grow = 1.2,
        unorm = u -> wnorm(u, norm_weights(cdata), Inf),
        rnorm = u -> wnorm(u, norm_weights(cdata), 1),
        solver_kwargs...
    )
    times = [i * period(voltages) for i in 0:nperiods]
    @info "Solving for $(voltages(0))V..."
    inival = solve(sys; inival = unknowns(esys), control = deepcopy(control), damp_initial = 0.1)
    result = CVSweepResult()
    allprogress = times[end] - times[begin]
    tprogress = 0
    @withprogress begin
        function pre(sol, t)
            working_electrode_voltage!(cdata,voltages(t))
        end

        function post(sol, oldsol, t, Δt)
            I = -VoronoiFVM.integrate(sys, sys.physics.breaction, sol; boundary = true)
            I_react = I[:, working_electrode(cdata)]
            I_we = -VoronoiFVM.integrate(sys, tf_we, sol, oldsol, Δt)
            I_bulk = -VoronoiFVM.integrate(sys, tf_bulk, sol, oldsol, Δt)
            push!(result.times, t)
            push!(result.voltages, voltages(t))
            push!(result.j_reaction, I_react)
            push!(result.j_we, I_we)
            push!(result.j_bulk, I_bulk)
            tprogress += abs(Δt)
            @logprogress tprogress / allprogress
        end

        function delta(sys, u, v, t, Δt)
            n = wnorm(u - v, norm_weights(cdata), Inf)
        end

        tsol = solve(
            sys;
            inival,
            times,
            control,
            pre,
            post,
            delta,
            store_all = store_solutions
        )
        if store_solutions
            result.tsol = tsol
        end
    end
    return result
end
