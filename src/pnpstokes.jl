"""
    $(TYPEDEF)

Stokes solver struct.
$(TYPEDFIELDS)
"""
struct PNPStokesSolver{Tflow, Tpnp}
    flowsolver::Tflow
    pnpsolver::Tpnp
end

"""
   flowsolver(grid; μ = 1, velospace = H1P2B)         
"""
function flowsolver end

function extended_unknowns end
function voltage! end
function node_pressure end
function node_velocity end
function fvm_velocities end
function flowplot end
function velocity_unknown end

"""
    PNPStokesSolver(;flowgrid, pnpgrid, μ, velospace, pnpdata, pnpbcond, pnpreaction, flowbcond, kwargs...)

Create Poisson-Nernst-Planck-Stokes-Solver.
"""
function PNPStokesSolver(;
        flowgrid = simplexgrid(0:0.5:1),
        pnpgrid = flowgrid,
        μ = 1,
        velospace = nothing,
        pnpdata = ElectrolyteData(),
        pnpbcond = (y, u, bnode, data) -> nothing,
        pnpreaction = (y, u, bnode, data) -> nothing,
        flowbcond = (flowsolver,) -> nothing,
        kwargs...
    )

    flowslv = flowsolver(flowgrid; μ, velospace)

    flowbcond(flowslv)

    pnpslv = PNPSystem(
        pnpgrid;
        celldata = pnpdata,
        bcondition = pnpbcond,
        reaction = pnpreaction,
        kwargs...
    )

    return PNPStokesSolver(flowslv, pnpslv)
end

function SciMLBase.solve(
        pnpssolver::PNPStokesSolver; damp0 = 0.1,
        embed = (data, λ) -> nothing,
        nembed = 0, niter = 20, damp_initial = 1, kwargs...
    )
    (; flowsolver, pnpsolver) = pnpssolver
    pnpstate = VoronoiFVM.SystemState(pnpsolver.vfvmsys)
    pnpdata = pnpsolver.vfvmsys.physics.data
    pnpgrid = pnpsolver.vfvmsys.grid
    (; iϕ, ip) = pnpdata

    flowsol = extended_unknowns(flowsolver)
    pnpsol = unknowns(pnpsolver.vfvmsys, inival = 0)
    t_pnp = 0.0
    t_stokes = 0.0
    t_project = 0.0

    q = zeros(num_nodes(pnpgrid))
    for iter in 1:niter
        if iter == 1
            inidamp = damp0
        else
            inidamp = damp_initial
        end
        if iter < nembed
            λ = iter / nembed
            embed(pnpdata, λ)
        else
            embed(pnpdata, 1.0)
        end
        t_pnp += @elapsed pnpsol = solve!(pnpstate; inival = pnpsol, damp_initial = inidamp, data = pnpdata, kwargs...)
        @views voltage!(flowsol, flowsolver, view(pnpsol, iϕ, :))
        chargedensity!(flowsol, flowsolver, chargedensity!(q, pnpsol, pnpdata))
        t_stokes += @elapsed solve!(flowsol, flowsolver; verbosity = -1)
        t_project += @elapsed evelo, bfvelo = fvm_velocities(flowsol, flowsolver; reconst = false)
        pnpdata.edgevelocity .= evelo
        pnpsol[ip,:].= node_pressure(flowsol, flowsolver)
        #        pnpdata.bfvelo.=bfvelo
    end
    @info "t_pnp=$(myround(t_pnp)), t_stokes=$(myround(t_stokes)) t_project=$(myround(t_project))"
    return pnpsol, flowsol
end
