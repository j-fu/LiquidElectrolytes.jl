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
function charge! end
function fvm_pressure end
function fvm_velocities end
function flowplot end
function velocity_unknown end

function PNPStokesSolver(;
                         flowgrid=simplexgrid(0:0.5:1),
                         pnpgrid=flowgrid,
                         μ=1,
                         velospace=nothing,
                         pnpdata=ElectrolyteData(),
                         pnpbcond=(y,u,bnode,data)->nothing,
                         pnpreaction=(y,u,bnode,data)->nothing,
                         flowbcond= (flowsolver,)-> nothing,
                         kwargs...)
    
    flowslv=flowsolver(flowgrid; μ, velospace)

    flowbcond(flowslv)

    pnpslv=PNPSystem(pnpgrid;
                     celldata=pnpdata,
                     bcondition = pnpbcond,
                     reaction = pnpreaction,
                     kwargs...)

    PNPStokesSolver(flowslv,pnpslv)
end

function SciMLBase.solve(pnpssolver::PNPStokesSolver; niter=10, kwargs...)
    (;flowsolver, pnpsolver)=pnpssolver
    pnpstate=VoronoiFVM.SystemState(pnpsolver)
    pnpdata=pnpsolver.physics.data
    pnpgrid=pnpsolver.grid
    (;iϕ)=pnpdata
   
    flowsol=extended_unknowns(flowsolver)
    pnpsol=unknowns(pnpsolver, inival=0)
    t_pnp=0.0
    t_stokes=0.0
    
    q=zeros(num_nodes(pnpgrid))
    for iter=1:niter
        t_pnp += @elapsed pnpsol=solve!(pnpstate; inival=pnpsol,  data=pnpdata, kwargs...)
        @views voltage!(flowsol, flowsolver,view(pnpsol,iϕ,:))
        charge!(flowsol,flowsolver,charge!(q,pnpsol, pnpdata))
        t_stokes +=@elapsed solve!(flowsol, flowsolver; verbosity=-1)
        evelo,bfvelo=fvm_velocities(flowsol, flowsolver; reconst =false)
        pnpdata.edgevelocity.=evelo
        pnpdata.pressure.= fvm_pressure(flowsol, flowsolver)
#        pnpdata.bfvelo.=bfvelo
    end
    @info "t_pnp=$(myround(t_pnp)), t_stokes=$(myround(t_stokes))"
    pnpsol, flowsol
end
