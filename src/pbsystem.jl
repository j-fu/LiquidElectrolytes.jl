"""
    pbspacecharge(φ, p, electrolyte)

Space charge expression for Poisson-Boltzmann
"""
function pbspacecharge(φ, p, electrolyte)
    c0_bulk, barc_bulk = c0_barc(electrolyte.c_bulk, electrolyte)
    pscaled = (p * electrolyte.pscale - electrolyte.p_bulk)
    c0 = c0_bulk * electrolyte.exp(-electrolyte.v0 * pscaled / (electrolyte.RT))
    sumyz = zero(p)
    sumyv = electrolyte.v0 * c0
    for α in 1:electrolyte.nc
        barv = electrolyte.v[α] + electrolyte.κ[α] * electrolyte.v0
        η_p = barv * pscaled
        η_φ = electrolyte.z[α] * electrolyte.F * (φ - electrolyte.ϕ_bulk)
        y = electrolyte.c_bulk[α] * rexp(-(η_φ + η_p) / (electrolyte.RT))
        sumyz += electrolyte.z[α] * y
        sumyv += barv * y
    end
    return electrolyte.F * sumyz / sumyv
end

"""
    pbreaction(f, u, node, electrolyte)

Reaction expression for Poisson-Boltzmann
"""
function pbreaction(f, u, node, electrolyte)
    iϕ, ip = 1, 2
    ## Charge density
    f[iϕ] = -pbspacecharge(u[iϕ], u[ip], electrolyte)
    f[ip] = 0
    return 
end

"""
    pbflux(f, u, edge, electrolyte)

Flux expression for Poisson-Boltzmann
"""
function pbflux(f, u, edge, electrolyte)
    iϕ, ip = 1, 2
    f[iφ] = electrolyte.ε * electrolyte.ε_0 * (u[iφ, 1] - u[iφ, 2])
    if iszero(electrolyte.v)
        qavg = 0
    else
        q1 = pbspacecharge(u[iφ, 1], u[ip, 1], electrolyte)
        q2 = pbspacecharge(u[iφ, 2], u[ip, 2], electrolyte)
        qavg = (q1 + q2) / 2
    end
    f[ip] = (u[ip, 1] - u[ip, 2]) + (u[iφ, 1] - u[iφ, 2]) * qavg / electrolyte.pscale 
    return 
end


"""
    PBSystem(grid;
             celldata=ElectrolyteData(),
             bcondition=default_bcondition,
             kwargs...)

Create VoronoiFVM system generalized Poisson-Boltzmann. Input:
- `grid`: discretization grid
- `celldata`: composite struct containing electrolyte data
- `bcondition`: boundary condition
- `kwargs`: Keyword arguments of VoronoiFVM.System
"""
function PBSystem(
        grid;
        celldata = ElectrolyteData(),
        bcondition = default_bcondition,
        kwargs...
    )

    return sys = VoronoiFVM.System(
        grid;
        data = celldata,
        flux = pbflux,
        reaction = pbreaction,
        bcondition,
        species = [1, 2],
        kwargs...
    )
end
