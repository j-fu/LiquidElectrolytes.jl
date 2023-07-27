"""
    splitz(range::AbstractRange)

If range contains zero, split it into two parts, one with values <=0 and one with values >=0.
Otherwise, return the range or its reverse, such that first value always is the one with the smallest absolute value.
"""
function splitz(range::AbstractRange)
    if range[1] >= 0
        return [range]
    elseif range[end] <= 0
        return [reverse(range)]
    else
        [0:(-step(range)):range[1], 0:step(range):range[end]]
    end
end

"""
    splitz(range::Vector)

Version of [`splitz(range::AbstractRange)`](@ref) for vectors.
"""
function splitz(range::Vector)
    if range[1] >= 0
        return [vcat([0.0], range)]
    elseif range[end] <= 0
        return [vcat([0.0], reverse(range))]
    else
        for i = 1:length(range)
            if range[i] ≈ 0.0
                return [reverse(range[1:i]), range[i:end]]
            elseif i > 1 && range[i - 1] < 0.0 && range[i] > 0.0
                return [vcat([0.0], reverse(range[1:(i - 1)])), vcat([0.0], range[i:end])]
            end
        end
    end
end

"""
    bulkbcondition(f,u,bnode,electrolyte)

Bulk boundary condition for electrolyte: set potential, pressure and concentrations to bulk values.
"""
function bulkbcondition(f, u, bnode, data; region = data.Γ_bulk)
    (; iϕ, ip, nc, ϕ_bulk, p_bulk, c_bulk) = data
    if bnode.region == region
        boundary_dirichlet!(f, u, bnode; species = iϕ, region, value = ϕ_bulk)
        boundary_dirichlet!(f, u, bnode; species = ip, region, value = p_bulk)
        for ic = 1:nc
            boundary_dirichlet!(f, u, bnode; species = ic, region, value = data.c_bulk[ic])
        end
    end
end

"""
    $(TYPEDEF)

Abstract simulation result.
"""
abstract type AbstractSimulationResult end

"""
    $(TYPEDEF)

Result data type for [`dlcapsweep`](@ref)

$(TYPEDFIELDS)
"""
struct DLCapSweepResult{Tv, Tc, Ts} <: AbstractSimulationResult
    """
    Vector of voltages
    """
    voltages::Tv

    """
    Vector of double layer capacitances
    """
    dlcaps::Tc

    """
    Vector of solutions
    """
    solutions::Ts
end


"""
    voltages_solutions(result)

Return a [`TransientSolution`](https://j-fu.github.io/VoronoiFVM.jl/stable/solutions/#Transient-solution) `tsol`
containing voltages (in `tsol.t`) and the corresponding stationary solutions (in `tsol.u`).
"""
function voltages_solutions end

voltages_solutions(r::DLCapSweepResult) = TransientSolution(r.solutions, r.voltages)

"""
    voltages_dlcaps(result)

Double layer capacitance curve as [`DiffEqArray`](https://docs.sciml.ai/RecursiveArrayTools/stable/array_types/#RecursiveArrayTools.DiffEqArray)
"""
voltages_dlcaps(r::DLCapSweepResult) = RecursiveArrayTools.DiffEqArray(r.dlcaps, r.voltages)

"""
           dlcapsweep(sys;voltages=(-1:0.1:1)*ufac"V",
                                  δ=1.0e-4,
                                  molarity=0.1*ufac"mol/dm^3",
                                  store_solutions=false,
                                  solver_kwargs...)

Calculate double layer capacitances for voltages given in `voltages`.
Returns a [`DLCapSweepResult`](@ref)

Assumptions:
- Only one double layer in the system - close to working electrode
- 1D domain
"""
function dlcapsweep(sys;
                    data=electrolytedata(sys),
                    inival = nothing,
                    iϕ = data.iϕ,
                    voltages = (-1:0.1:1) * ufac"V",
                    δ = 1.0e-4,
                    molarity = 0.1 * ufac"mol/dm^3",
                    store_solutions = false,
                    solver_kwargs...)
    ranges = splitz(voltages)
    vplus = zeros(0)
    cdlplus = zeros(0)
    vminus = zeros(0)
    cdlminus = zeros(0)
    splus = []
    sminus = []

    data.ϕ_we = 0

    data.c_bulk .= molarity

    if isnothing(inival)
        inival=pnpunknowns(sys)
    end 

    inival = solve(sys; inival, verbose="n",damp_initial=0.1)
    allprogress = sum(length, ranges)
    ϕprogress = 0

    function show_error(u, δ)
        @show u[iϕ, 1:5]
        @show u[ip, 1:5]
        @error "bailing out at δ=$(δ) ϕ_we=$(data.ϕ_we)V, molarity=$(molarity)"
    end

    control = VoronoiFVM.SolverControl(; max_round = 3, tol_round = 1.0e-9, solver_kwargs...)
    @withprogress for range in ranges
        sol = inival
        success = true
        for ϕ in range
            try
                data.ϕ_we = ϕ
                sol = solve(sys; inival = sol, control)
            catch e
                println(e)
                show_error(sol, 0)
                success = false
            end

            sol0 = sol

            if !success
                break
            end
            Q = integrate(sys, sys.physics.reaction, sol)

            try
                data.ϕ_we = ϕ + δ
                sol = solve(sys; inival = sol, control)
            catch e
                println(e)
                show_error(sol, δ)
                success = false
            end
            if !success
                break
            end
            Qδ = integrate(sys, sys.physics.reaction, sol)
            cdl = (Qδ[iϕ] - Q[iϕ]) / δ

            if range[end] > range[1]
                store_solutions ? push!(splus, sol) : nothing
                push!(vplus, ϕ)
                push!(cdlplus, cdl)
            else
                store_solutions ? push!(sminus, sol) : nothing
                push!(vminus, ϕ)
                push!(cdlminus, cdl)
            end
            ϕprogress += 1
            @logprogress ϕprogress / allprogress
        end
    end

    volts = vcat(reverse(vminus), vplus)
    cdls = vcat(reverse(cdlminus), cdlplus)
    GC.gc()
    DLCapSweepResult(volts, cdls, store_solutions ? vcat(reverse(sminus), splus) : nothing)
end

"""
    $(TYPEDEF)

Result data type for [`ivsweep`](@ref)

$(TYPEDFIELDS)
"""
struct IVSweepResult{Tv, Twe, Tbulk, Ts} <: AbstractSimulationResult
    """
    Vector of voltages
    """
    voltages::Tv

    """
    Working electrode molar reaction rates
    """
    j_we::Twe

    """
    Bulk molar fluxes
    """
    j_bulk::Tbulk

    """
    Vector of solutions
    """
    solutions::Ts
end

voltages_solutions(r::IVSweepResult) = TransientSolution(r.solutions, r.voltages)

"""
    voltages_currents(result,ispec)

Voltage- working electrode current curve for species as [`DiffEqArray`](https://docs.sciml.ai/RecursiveArrayTools/stable/array_types/#RecursiveArrayTools.DiffEqArray)
"""
voltages_currents(r::IVSweepResult, ispec) = RecursiveArrayTools.DiffEqArray(currents(r, ispec), r.voltages)

"""
    currents(result,ispec)

Working electrode current  for species `ispec`.
"""
currents(r::IVSweepResult, ispec) = [ph"F" * j[ispec] for j in r.j_we]

"""
     ivsweep(
          sys;
          voltages = (-0.5:0.1:0.5) * ufac"V",
          store_solutions = false,
          solver_kwargs...,
          )

Calculate molar reaction rates and bulk flux rates for each voltage in `voltages`.
"""
function ivsweep(sys;
                 voltages = (-0.5:0.1:0.5) * ufac"V",
                 store_solutions = false,
                 solver_kwargs...)
    ranges = splitz(voltages)
    F = ph"N_A*e"
    data = sys.physics.data

    factory = VoronoiFVM.TestFunctionFactory(sys)
    tf_bulk = testfunction(factory, [data.Γ_we], [data.Γ_bulk])

    iplus = []
    iminus = []
    fplus = []
    fminus = []
    vminus = zeros(0)
    vplus = zeros(0)
    sminus = []
    splus = []
    data = electrolytedata(sys)
    data.ϕ_we = 0
    control = SolverControl(;
                            verbose = true,
                            handle_exceptions = true,
                            Δp_min = 1.0e-3,
                            Δp = 1.0e-2,
                            Δp_grow = 1.2,
                            Δu_opt = 1.0e-2,
                            unorm = u -> wnorm(u, data.weights, Inf),
                            rnorm = u -> wnorm(u, data.weights, 1),
                            solver_kwargs...)
    iϕ = data.iϕ
    @info "Solving for 0V..."
    inival = solve(sys; inival = pnpunknowns(sys), control)

    allprogress = voltages[end] - voltages[begin]
    ϕprogress = 0
    for range in ranges
        @info "IV sweep from $(range[1])V to $(range[end])V..."
        dir = range[end] > range[1] ? 1 : -1

        psol = nothing
        @withprogress begin
            function pre(sol, ϕ)
                data.ϕ_we = dir * ϕ
            end

            function post(sol, oldsol, ϕ, Δϕ)
                I_react = -integrate(sys, sys.physics.breaction, sol; boundary = true)[:,
                                                                                       data.Γ_we]
                I_bulk = -integrate(sys, tf_bulk, sol)
                if dir > 0
                    push!(vplus, data.ϕ_we)
                    push!(iplus, I_react)
                    push!(fplus, I_bulk)
                else
                    push!(vminus, data.ϕ_we)
                    push!(iminus, I_react)
                    push!(fminus, I_bulk)
                end
                ϕprogress += abs(Δϕ)
                @logprogress ϕprogress / allprogress
            end

            function delta(sys, u, v, t, Δt)
                n = wnorm(u, data.weights, Inf) * data.v0
            end

            psol = solve(sys;
                         inival,
                         embed = dir * range,
                         control,
                         pre,
                         post,
                         delta,
                         store_all = store_solutions)
        end

        if dir == 1
            splus = psol.u
        else
            sminus = psol.u

            popfirst!(iminus)
            popfirst!(vminus)
            popfirst!(fminus)
            if store_solutions
                popfirst!(sminus)
            end
        end
    end

    IVSweepResult(vcat(reverse(vminus), vplus),
                  vcat(reverse(iminus), iplus),
                  vcat(reverse(fminus), fplus),
                  store_solutions ? vcat(reverse(sminus), splus) : nothing)
end
