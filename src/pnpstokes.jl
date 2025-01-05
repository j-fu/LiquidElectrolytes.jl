struct PNPStokesSolver
    flowsolver::FlowSolver
    pnpsolver::VoronoiFVM.System
end

function PNPStokesSolver(;
                         flowgrid=simplexgrid(0:0.5:1),
                         pnpgrid=flowgrid,
                         μ=1,
                         velospace = H1P2B,
                         pnpdata=ElectrolyteData(),
                         pnpbcond=(y,u,bnode,data)->nothing,
                         pnpreaction=(y,u,bnode,data)->nothing,
                         flowbcond= (flowsolver,)-> nothing,
                         kwargs...)
    
    flowslv=FlowSolver(flowgrid; μ, velospace)
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
    (;iφ)=pnpdata
   
    flowsol=extended_unknowns(flowsolver)
    pnpsol=unknowns(pnpsolver, inival=0)
    
    q=zeros(num_nodes(pnpgrid))
    for iter=1:niter
        @info "iter=$(iter)"
        pnpsol=solve!(pnpstate; inival=pnpsol,  data=pnpdata, kwargs...)
        @views voltage!(flowsol, flowsolver,view(pnpsol,pnpdata.iϕ,:))
        charge!(flowsol,flowsolver,charge!(q,pnpsol, pnpdata))
        solve!(flowsol, flowsolver)
        evelo,bfvelo=fvm_velocities(flowsol, flowsolver; reconst =false)
        pnpdata.evelo.=evelo
#        pnpdata.bfvelo.=bfvelo
    end
    pnpsol, flowsol
end
