module LiquidElectrolytesExtendableFEMExt

using SciMLBase
using ExtendableFEM
using ExtendableGrids
using LiquidElectrolytes
using VoronoiFVM

    
function kernel_stokes_cartesian!(result, u_ops, qpinfo)
    μ = qpinfo.params[1]
    ∇u, p = view(u_ops, 1:4), view(u_ops, 5)
    result[1] = μ * ∇u[1] - p[1]
    result[2] = μ * ∇u[2]
    result[3] = μ * ∇u[3]
    result[4] = μ * ∇u[4] - p[1]
    result[5] = -(∇u[1] + ∇u[4])
    return nothing
end

function kernel_stokes_cylindrical!(result, u_ops, qpinfo)
    μ = qpinfo.params[1]
    r = qpinfo.x[1]
    u, ∇u, p = view(u_ops, 1:2), view(u_ops, 3:6), view(u_ops, 7)
    result[1] = μ / r * u[1] - p[1]
    result[2] = 0
    result[3] = μ * r * ∇u[1] - r * p[1]
    result[4] = μ * r * ∇u[2]
    result[5] = μ * r * ∇u[3]
    result[6] = μ * r * ∇u[4] - r * p[1]
    result[7] = -(r * (∇u[1] + ∇u[4]) + u[1])
    return nothing
end

function kernel_lorentz!(result, args, qpinfo)
    params = qpinfo.params
    q = view(args, 3)
    result .= view(args, 1:2) * q[1]
    return nothing
end


struct FlowSolver
    grid::ExtendableGrid
    problem::ProblemDescription
    u::Unknown
    p::Unknown
    φ::Unknown
    q::Unknown
    FE_u::FESpace
    FE_p::FESpace
    FE_ϕ::FESpace
    FE_q::FESpace
    T_fereconst
    iu::Int
    ip::Int
    iφ::Int
    iq::Int

end

function LiquidElectrolytes.flowsolver(grid::ExtendableGrid; μ = 1, velospace = H1P2B)
    u = Unknown("u"; name = "velocity")
    p = Unknown("p"; name = "pressure")
    φ = Unknown("ϕ"; name = "voltage")
    q = Unknown("q"; name = "charge")
    
    problem = ProblemDescription("incompressible Stokes problem")
    assign_unknown!(problem, u)
    assign_unknown!(problem, p)
    
    stokes_op = BilinearOperator(
        kernel_stokes_cartesian!,
        [grad(u), id(p)], [grad(u), id(p)],
        bonus_quadorder = 2, store = false,
        params = [μ]
    )
    
    lorentz_op = LinearOperator(
        kernel_lorentz!,
        # [apply(1, Reconstruct{HDIVRT0{2}, Identity})],
        [id(1)],
        [grad(3), id(4)];
        bonus_quadorder = 3
    )
    
    assign_operator!(problem, stokes_op)
    assign_operator!(problem, lorentz_op)
    if velospace == H1P2B
        FE_u = FESpace{H1P2B{2, 2}}(grid)
        FE_p = FESpace{L2P1{1}}(grid, broken = true)
        T_fereconst = HDIVBDM2{2}
    elseif velospace == H1BR
        FE_u = FESpace{H1BR{2}}(grid)
        FE_p = FESpace{L2P0{1}}(grid, broken = true)
        T_fereconst = HDIVRT0{2}
    else
        error("Wrong velospace, choose between H1P2B and H1BR")
    end
    FE_φ = FESpace{H1P1{1}}(grid)
    FE_q = FESpace{H1P1{1}}(grid)
    
    return FlowSolver(
        grid,
        problem,
        u, p, φ, q,
        FE_u, FE_p, FE_φ, FE_q,
        T_fereconst,
        1, 2, 3, 4
    )
end



ExtendableFEM.assign_operator!(fs::FlowSolver, op) = assign_operator!(fs.problem, op)

function  LiquidElectrolytes.extended_unknowns(fs::FlowSolver)
    return FEVector(
        [
            fs.FE_u,
            fs.FE_p,
            fs.FE_ϕ,
            fs.FE_q,
        ],
        tags = [fs.u, fs.p, fs.φ, fs.q]
    )
end

LiquidElectrolytes.velocity_unknown(fs::FlowSolver) = fs.u

function  LiquidElectrolytes.voltage!(sol::FEVector, fs::FlowSolver, voltage)
    return view(sol[fs.iφ]) .= voltage
end
function  LiquidElectrolytes.charge!(sol::FEVector, fs::FlowSolver, charge)
    return view(sol[fs.iq]) .= charge
end


function SciMLBase.solve!(sol::FEVector, fs::FlowSolver; kwargs...)
    return ExtendableFEM.solve(fs.problem, [fs.FE_u, fs.FE_p]; init = sol, kwargs...)
end

function SciMLBase.solve(
        fs::FlowSolver;
        voltage = zeros(num_nodes(fs.grid)),
        charge = zeros(num_nodes(fs.grid)),
        kwargs...
    )
    sol = extended_unknowns(fs)
    voltage!(sol, fs, voltage)
    charge!(sol, fs, charge)
    return solve!(sol, fs; kwargs...)
end

function LiquidElectrolytes.fvm_velocities(
        sol::FEVector, fs::FlowSolver;
        fvmgrid = fs.grid,
        reconst = true,
        interpolation_eps = 1.0e-9
    )
    flowgrid = fs.grid
    FES_reconst = FESpace{fs.T_fereconst}(flowgrid)
    R = FEVector(FES_reconst)
    cylindrical_grid = false
    if reconst
        if cylindrical_grid
            lazy_interpolate!(
                R[1], sol, [id(fs.u)]; postprocess = multiply_r,
                bonus_quadorder = 2, eps = interpolation_eps
            )
            femvelo = R[1]
        else
            lazy_interpolate!(
                R[1], sol, [id(fs.u)];
                bonus_quadorder = 2, eps = interpolation_eps
            )
            femvelo = R[1]
        end
    else
        femvelo = sol[1]
    end
    evelo = VoronoiFVM.edgevelocities(fvmgrid, femvelo; reconst)
    bfvelo = VoronoiFVM.bfacevelocities(fvmgrid, femvelo; reconst)
    return evelo, bfvelo
end

LiquidElectrolytes.fvm_pressure(flowsol, flowsolver) = view(nodevalues(flowsol[flowsolver.p]), 1,:)

function LiquidElectrolytes.flowplot(sol::FEVector, fs::FlowSolver; kwargs...)
    return ExtendableFEM.plot([id(fs.u), id(fs.p)], sol; kwargs...) |> reveal
end



#############################################
# Deprecated

function stokes_operator(u, p; cyl = false, μ = 1)
    if cyl
        return BilinearOperator(
            kernel_stokes_cylindrical!,
            [id(u), grad(u), id(p)],
            bonus_quadorder = 2, store = false,
            params = [μ]
        )
    else
        return BilinearOperator(
            kernel_stokes_cartesian!,
            [grad(u), id(p)],
            bonus_quadorder = 2, store = false,
            params = [μ]
        )
    end
end

function stokes_problem(u, p; cyl = false, μ = 1)
    op = stokes_operator(u, p; cyl, μ)
    problem = ProblemDescription("incompressible Stokes problem")
    assign_unknown!(problem, u)
    assign_unknown!(problem, p)
    assign_operator!(problem, op)
    return problem
end

function stokes_space(grid)
    u = Unknown("u"; name = "velocity")
    p = Unknown("p"; name = "pressure")
    FE_u = H1BR{2}
    FE_p = L2P0{1}
    return u, p, [FESpace{FE_u}(grid), FESpace{FE_p}(grid; broken = true)]
end


function multiply_r!(result, input, qpinfo)
    x = qpinfo.x
    result .= input * x[1]
    return nothing
end

end
